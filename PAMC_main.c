#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "dSFMT.h"
#include "./Lorenz96/model_L96.h"

int get_Lidx(float Lfrac, int *Lidx);
void get_Lidx_from_file(int L, int *Lidx, char *fname);
void loadSpecs(double *xob, double *xold, int *Lidx, double *p);
void loadMeasurements(int L, int Nob, double *xob, char *fname);
void loadInitialPath(double *xold, char *fname);
double calc_measErr(double *y, double *x);
double calc_modelErr(double *x, double *p);
double calc_A(double *y, double *x, double *p);
int calc_dA_F(int np, double xnew_md, double *x, double *p);
double calc_pointMeasErr(int m, int d, double *y, double *x);
double calc_pointModelErr(int m, int d, double *x, double *p);
double calc_dA(int m, int d, double xnew_md, double *y, double *x, double *p);
double RK2(int m,int d,double *x_m,double *p);
void writeToFiles(FILE *fpA, FILE *fpP, double *y, double *x, double *p);
void writeFinalPath(FILE *fpX, double *x);


double Rm[D];
double Rf[D];
int Lidx[D];
double delta_x, xacc_expected;
char *PATH;
int IniSeed;


int main(int argc, char* argv[]){
    if(argc!=3){
        printf("Usage: ./PAMC_main.c No.ofInitialConditions <path>\n");
        exit(1);
    }
    PATH = argv[2];
    IniSeed = atoi(argv[1]);//Number for different initial conditions
    printf("Starting initial condition %d",IniSeed);
    double xob[M*D],xold[M*D],xblock[M*D],xave[M*D];
    double pold[NP],pblock[NP],pave[NP];
    double xnew, dA0, Action, Pacc;
    double SEED = 3154+IniSeed;
    int m,d,i,np,beta;

    //Variables used for updating step size in the initialization phase

    double delta[M*D];//X step size
    double delta_F[NP];//P step size
    double updateSize = 100;
    int Nacc[M*D];
    int Nacc_F[NP];

    dsfmt_t dsfmt;
    dsfmt_init_gen_rand(&dsfmt, SEED);


    for(np=0;np<NP;np++){
        delta_F[np] = delta_0;
        Nacc_F[np] = 0;
    }

    printf("Model Parameter:M = %d, D = %d, NP = %d\n", M, D, NP);

    // Loading specs from the file specs.txt
    loadSpecs(xob, xold, Lidx, pold);

    // Initialization of some parameters and state variables
    for(m=0;m<M;m++){
        for(d=0;d<D;d++){
            i = m*D+d;
            delta[i] = delta_0;
            Nacc[i] = 0;
            xblock[i] = 0;
            xave[i] = xold[i];
        }
    }

    for(np=0;np<NP;np++){
        delta_F[np] = delta_0;
        Nacc_F[np] = 0;
        pblock[np] = 0;
        pave[np] = pold[np];
    }

    int Nit = N_ini + Bsize*Nblock;
    int iIt, count=0;

    // Open files to write results. action, final path, parameter est
    char fn1[255],fn2[255],fn3[255];
    sprintf(fn1, "%s/outputs/action_%d.dat", PATH, IniSeed);
    sprintf(fn2, "%s/outputs/path_%d.dat", PATH, IniSeed);
    sprintf(fn3, "%s/outputs/parameter_%d.dat", PATH, IniSeed);
    FILE *fp1 = fopen(fn1,"w");
    FILE *fp2 = fopen(fn2,"w");
    FILE *fp3 = fopen(fn3,"w");

    for(beta=0;beta<beta_max;beta++){
        // Initialization needed after each beta
        for(m=0;m<M;m++){
            for(d=0;d<D;d++){
                i = m*D+d;
                xold[i] = xave[i];
                Nacc[i] = 0;
                xblock[i] = 0;
                xave[i] = 0;
            }
        }
        for(np=0;np<NP;np++){
            pold[np] = pave[np];
            Nacc_F[np] = 0;
            pblock[np] = 0;
            pave[np] = 0;
        }
        Action = calc_A(xob, xold, pold);
        printf("Action: %f\n", Action);
        // Perturb all variables including paramters Nit times
        for(iIt=0;iIt<Nit;iIt++){
            // Perturbing state variables
            for(m=0;m<M;m++){
                for(d=0;d<D;d++){
                    i = m*D+d;
                    xnew = xold[i] + delta[i]*(2.0*dsfmt_genrand_open_open(&dsfmt)-1.0);
                    dA0 = calc_dA(m,d,xnew,xob,xold,pold);
                    if(dA0<0){
                        Nacc[i]++;
                        xold[i] = xnew;
                        Action += dA0;
                    }
                    else{
                        Pacc = exp(-dA0);
                        if(dsfmt_genrand_open_open(&dsfmt) <= Pacc){
                            Nacc[i]++;
                            xold[i] = xnew;
                            Action += dA0;
                        }//end if
                    }//end else
                }//end for d
            }//end for m

            // Perturbing parameters
            for(np=0;np<NP;np++){
                xnew = pold[np] + delta_F[np]*(2.0*dsfmt_genrand_open_open(&dsfmt)-1.0);
                dA0 = calc_dA_F(np, xnew, xold, pold);
                if(dA0<0){
                    Nacc_F[np]++;
                    pold[np] = xnew;
                    Action += dA0;
                }
                else{
                    Pacc = exp(-dA0);
                    if(dsfmt_genrand_open_open(&dsfmt) <= Pacc){
                        Nacc_F[np]++;
                        pold[np] = xnew;
                        Action += dA0;
                    }//end if
                }//end else
            }//end for np

            //Initialization phase
            if(iIt < N_ini){
                if((iIt+1)%((int)updateSize)==0){
                    for(m=0;m<M;m++){
                        for(d=0;d<D;d++){
                            i = m*D+d;
                            delta[i] = delta[i]*(1.0 + delta_x*(Nacc[i]/updateSize-xacc_expected));
                            Nacc[i] = 0;
                        }//end for d
                    }//end for m
                    for(np=0;np<NP;np++){
                        delta_F[np] = delta_F[np]*(1.0 + delta_x*(Nacc_F[np]/updateSize-xacc_expected));
                        Nacc_F[np] = 0;
                    }//end for np
                }//end if
            }//end if initialization

            if(iIt == (N_ini-1)){
                printf("Initialization phase done.\n");
                printf("delta[0] = %f, delta_F[0] = %f\n",delta[0], delta_F[0]);
            }//end if iIt==N_ini

            // Main MC steps
            if(iIt>=N_ini){
                for(m=0;m<M;m++){
                    for(d=0;d<D;d++){
                        i = m*D+d;
                        xblock[i] += xold[i]/Bsize;
                    }//end for d
                }//end for m
                for(np=0;np<NP;np++){
                    pblock[np] += pold[np]/Bsize;
                }//end for np
                count++;

                // Add up the average of small blocks
                if(count%Bsize==0){
                    for(m=0;m<M;m++){
                        for(d=0;d<D;d++){
                            i = m*D+d;
                            xave[i] += xblock[i]/Nblock;
                            xblock[i] = 0;
                        }//end for d
                    }//end for m
                    for(np=0;np<NP;np++){
                        pave[np] += pblock[np]/Nblock;
                        pblock[np] = 0;
                    }//end for np
                }//end if count%Bsize==0
            }//end if beta==beta_max
        }//end for iIt

        printf("Count: %d\n",count);

        printf("Writing beta=%d to files.\n", beta);
        writeToFiles(fp1, fp3, xob, xave, pave);
        if(beta==(beta_max-1)){
            writeFinalPath(fp2, xave);
        }

        printf("beta=%d of %d done.\n",beta,beta_max);
        for(d=0;d<D;d++){
            Rf[d] = Rf[d]*alpha;
        }
    }//end for beta
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);

    return 0;
}//end main

//forward integrator, 2nd order Runge-Kutta
double RK2(int m, int d, double *x_m, double *p){
    double a1[D],x_out[D];
    int d_temp,step;
    if(nstep==1){
        for(d_temp=0;d_temp<D;d_temp++){
            a1[d_temp] = x_m[d_temp] + 0.5*dt*f(d_temp,m,x_m,p);
        }//end for d_temp
        return x_m[d] + dt*f(d,m,a1,p);
    }//end if
    else{
        for(d_temp=0;d_temp<D;d_temp++){
            x_out[d] = x_m[d];
        }//end for d_temp
        for(step=0;step<nstep;step++){
            for(d_temp=0;d_temp<D;d_temp++){
                a1[d_temp] = x_out[d_temp] + 0.5*dt*f(d_temp,m,x_out,p);
            }//end for d_temp
            for(d_temp=0;d_temp<D;d_temp++){
                x_out[d_temp] += dt*f(d_temp,m,a1,p);
            }//end for d_temp
        }//end for step
        return x_out[d];
    }//end else
}

//calculate the measurement error term in action
double calc_measErr(double *y, double *x){
    double error = 0;
    int d, m, i;
    for(m=0;m<M;m++){
        for(d=0;d<D;d++){
            if(Lidx[d]){
                i = m*D+d;
                error += 0.5*Rm[d]*(y[i]-x[i])*(y[i]-x[i]);
                //printf("%d,%d,%d,%f,%f,%f\n",m,d,Lidx[d],y[i],x[i],error);
            }
        }
    }
    return error;
}

//calculate the model error term in action
double calc_modelErr(double *x, double *p){
    double fx, error = 0;
    int d, m, i;
    for(m=1;m<M;m++){
        for(d=0;d<D;d++){
            i = m*D+d;
            fx = RK2(m-1,d,&x[m*D-D],p);
            error += 0.5*Rf[d]*(x[i]-fx)*(x[i]-fx);
            //printf("%f,%d,%d\n",fx,m,d);
        }
    }
    return error;
}


//calculate total action, simply summing up the previous 2 functions
double calc_A(double *y, double *x, double *p){
    double meE;
    double moE;
    meE = calc_measErr(y, x);
    moE = calc_modelErr(x, p);
    //printf("Measurement Error: %f\n",meE);
    //printf("Model Error: %f\n",moE);
    return meE+moE;
}

//calculate change in measurement error term due to perturbation in MC step
double calc_pointMeasErr(int m, int d, double *y, double *x){
    int i = m*D+d;
    if(Lidx[d]){
        return 0.5*Rm[d]*(y[i]-x[i])*(y[i]-x[i]);
    }
    else{
        return 0;
    }
}

//calculate change in model error term due to perturbation in MC step
double calc_pointModelErr(int m, int d, double *x, double *p){
    double fx;
    int i = m*D+d;
    fx = RK2(m,d,&x[m*D],p);
    return 0.5*Rf[d]*(x[i+D]-fx)*(x[i+D]-fx);
}

//calculate change in total action due to perturbation in MC step
double calc_dA(int m, int d, double xnew_md, double *y, double *x, double *p){
    double dA[2] = {0,0}; //dA[0] is the old pointwise action, dA[1] is the new one.
    double temp;
    int dd,i;
    i = m*D+d;
    int itemp, dtemp;
    double xm[D];
    dA[0] = calc_pointMeasErr(m,d,y,x);
    //printf("measEold: %lf\t",dA[0]);
    //perturb x_m,d, calculate change in x_m,d - f(x_m-1) if m>0
    if(m>0){
        dA[0] += calc_pointModelErr(m-1,d,x,p);
    }
    //perturb x_m,d, calculate change in sum{x_m+1,i - f(x_m)}^2 if m<M
    if(m<M-1){
        for(dd=0;dd<D;dd++){
            dA[0] += calc_pointModelErr(m,dd,x,p);
        }
    }
    temp = x[i];
    x[i] = xnew_md;
    dA[1] = calc_pointMeasErr(m,d,y,x);
    //printf("measEnew: %lf\n",dA[1]);
    if(m>0){
        dA[1] += calc_pointModelErr(m-1,d,x,p);
    }
    if(m<M-1){
        for(dd=0;dd<D;dd++){
            dA[1] += calc_pointModelErr(m,dd,x,p);
        }
    }
    x[i] = temp;
    return dA[1]-dA[0];
}

//calculate change in action due to perturbation in parameter
int calc_dA_F(int np, double xnew_md, double *x, double *p){
    double dA[2] = {0,0}; //dA[0] is the old pointwise action, dA[1] is the new one.
    double temp;
    dA[0] = calc_modelErr(x, p);
    temp = p[np];
    p[np] = xnew_md;
    dA[1] = calc_modelErr(x, p);
    p[np] = temp;
    return dA[1]-dA[0];
}

//load specs from specs.txt
void loadSpecs(double *xob, double *xold, int *Lidx, double *p){
    int d, Nob, out, i1, i2, i3, i;
    int L = 0;
    float f1, f2, f3, Lfrac;
    char key1, key2;
    char line[127], fname[127], fname2[127], name[127], RmfFile[127];
    FILE *sp;
    FILE *Rm_Rf;
    for (d=0;d<D;d++){
        Rm[d] = 0;
        Rf[d] = 0;
    }
    sprintf(name, "%s/specs.txt", PATH);
    sp = fopen(name,"r");
    if(sp==NULL){
        printf("File %s is missing\n",name);
        exit(1);
    }
    printf("Reading file %s\n", name);
    while(fgets(line,127,sp)!=NULL){
        if(line[0]!='#' && strlen(line)>2){
            key1 = line[0];
            key2 = line[1];

            // Annealing parameter, Rm and Rf0. Same for all variables
            if(key1=='A' && key2=='P'){
                out = sscanf(&line[2], "%f %f", &f1,&f2);
                for(d=0;d<D;d++){
                    Rm[d] = f1;
                    Rf[d] = f2;
                }
            }

            // Annealing parameter from file, in case that they are different
            // for different variables.
            if(key1=='A' && key2=='F'){
                out = sscanf(&line[2],"%s %s", fname,fname2);

                sprintf(RmfFile,"%s/%s", PATH,fname);
                printf("Loading Rm from %s.\n",RmfFile);
                Rm_Rf = fopen(RmfFile,"r");
                for(i=0;i<D;i++){
                    out = fscanf(Rm_Rf,"%f",&f1);
                    Rm[i] = f1;
                }
                fclose(Rm_Rf);

                sprintf(RmfFile,"%s/%s", PATH,fname2);
                printf("Loading Rf from %s.\n", RmfFile);
                Rm_Rf = fopen(RmfFile,"r");
                for(i=0;i<D;i++){
                    out = fscanf(Rm_Rf,"%f",&f1);
                    Rf[i] = f1;
                }
                fclose(Rm_Rf);
            }

            // Parameters for the step size adjustment during initialization phase
            if(key1=='S' && key2=='S'){
                out = sscanf(&line[2], "%f %f", &f1,&f2);
                delta_x = f1;
                xacc_expected = f2;
            }

            // Load measurement index. This option only given the fraction/percentage
            // of measured variables. Then assign measurements uniformly spaced.
            // Use this option for twin experiments or test case.
            if(key1=='L' && key2=='F'){
                if(L!=0){
                    printf("Observed index already set. Don't use both LF and OD keywords in spec.txt.\n");
                    exit(1);
                }
                out = sscanf(&line[2], "%f", &f1);
                Lfrac = f1;
                L = get_Lidx(Lfrac, Lidx);
            }

            // Load measurment index. This option will read measured index from the file.
            // Use this option for actual data set might be more convenient.
            if(key1=='O' && key2=='D'){
                if(L!=0){
                    printf("Observed index already set. Don't use both LF and OD keywords in spec.txt.\n");
                    exit(1);
                }
                out = sscanf(&line[2], "%d %s", &i1,fname);
                L = i1;
                get_Lidx_from_file(L,Lidx,fname);
            }

            // Load measurements. Nob can be equal to L (for real measurements)
            // or equal to D (for twin experiments or test cases)
            if(key1=='M' && key2=='P'){
                out = sscanf(&line[2], "%d %s", &i1,fname);
                Nob = i1;
                if(L==0){
                    printf("Observed index not set yet.\n");
                    printf("Set the observed index using LF or OD in spec.txt.\n");
                    exit(1);
                }
                loadMeasurements(L, Nob, xob, fname);
            }

            // Load inital paths from file.
            if(key1=='I' && key2=='P'){
                out = sscanf(&line[2], "%s", fname);
                loadInitialPath(xold, fname);
            }

            // Load inital guess on parameters. Right now option "PA" only allow
            // single parameter estimation. Will update later to load from file
            // for multiple parameter initial guesses.
            if(key1=='P' && key2=='A'){
                if((int)(sizeof(p)/sizeof(p[0]))!=1){
                    printf("More than 1 parameters. Use a file to load initial guesses for all poldmeters.\n");
                    exit(1);
                }
                out = sscanf(&line[2], "%f", &f1);
                p[0] = f1;
            }
        }
    }
    fclose(sp);
}

int get_Lidx(float Lfrac, int *Lidx){
    if(Lfrac>1){
        printf("Invalid Lfrac.\n");
        exit(1);
    }
    float inverse_Lfrac = 1.0/Lfrac;
    float temp;
    int indx;
    int Ltemp = 0;
    for(temp=0;temp<D;temp = temp+inverse_Lfrac){
        indx = (int)temp;
        Lidx[indx] = 1;
        Ltemp  = Ltemp +1;
    }
    return Ltemp;
}

void get_Lidx_from_file(int L, int *Lidx, char *fname){
    int l,temp,out;
    FILE *fLidx;
    char name[127];
    sprintf(name, "%s/%s", PATH,fname);
    fLidx = fopen(name,"r");
    if(fLidx==NULL){
        printf("File %s for observed variable index is missing\n",name);
        exit(1);
    }
    for(l=0;l<L;l++){
        out = fscanf(fLidx," %d",&temp);
        Lidx[temp] = 1;
    }
    fclose(fLidx);
}

void loadMeasurements(int L, int Nob, double *xob, char *fname){
    int l,m,i,out;
    double f1;
    FILE *fob;
    char name[127];
    if(Nob==D){
        sprintf(name, "%s/%s", PATH,fname);
        fob = fopen(name,"r");
        if(fob==NULL){
            printf("File %s is missing\n",name);
            exit(1);
        }
        printf("Loading observations from %s\n", name);
        for(m=0;m<M;m++){
            for(l=0;l<Nob;l++){
                i = m*D+l;
                out = fscanf(fob," %lf",&f1);
                xob[i] = f1;
            }
        }
    }
    else if(Nob==L){
        sprintf(name, "%s/%s", PATH,fname);
        fob = fopen(name,"r");
        if(fob==NULL){
            printf("File %s is missing\n",name);
            exit(1);
        }
        printf("Loading observations from %s\n", name);
        for(m=0;m<M;m++){
            for(l=0;l<Nob;l++){
                if(Lidx[l]){
                    i = m*D+l;
                    out = fscanf(fob," %lf",&f1);
                    xob[i] = f1;
                }
            }
        }
    }
    fclose(fob);
}

void loadInitialPath(double *xold, char *fname){
    int m,d,i,out;
    double f1;
    FILE *fxold;
    char name[127];
    sprintf(name, "%s/%s_%d.dat", PATH,fname,IniSeed);
    fxold = fopen(name,"r");
    if(fxold==NULL){
        printf("File %s for initial path is missing\n",name);
        exit(1);
    }
    printf("Loading initial path from %s\n", name);
    for(m=0;m<M;m++){
        for(d=0;d<D;d++){
            i = m*D+d;
            out = fscanf(fxold,"%lf",&f1);
            xold[i] = f1;
        }
    }
    fclose(fxold);
}

// Write action and parameter estimation to file at the end of each beta.
void writeToFiles(FILE *fpA, FILE *fpP, double *y, double *x, double *p){
    int np;
    double modelE, measE;
    modelE = calc_modelErr(x,p);
    measE = calc_measErr(y,x);
    fprintf(fpA, "%e\t%e\t%e\n", (modelE+measE),measE,modelE);
    for(np=0;np<NP;np++){
        fprintf(fpP, "%f\t", p[np]);
        fprintf(fpP, "\n");
    }
}

// Write the path of variables to file at the last beta.
void writeFinalPath(FILE *fpX, double *x){
    int m,d,i;
    for(m=0;m<M;m++){
        for(d=0;d<D;d++){
            i = m*D+d;
            fprintf(fpX, "%f\t", x[i]);
        }
        fprintf(fpX, "\n");
    }
}
