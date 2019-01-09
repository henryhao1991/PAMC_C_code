#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#include "model_swe.h"

/* This is the file to run shallow water equations forward in time, in order to
test if the model is coded right or not. */

void RK4(int m, double *x_m, double *x_mp1, double *p, int m_step);
void loadInitialPath(double *xold, char *fname);

int main(int argc, char* argv[]){

  double xold[D], x_mp1[D];
  int m,d;
  int m_step = 10;
  double p[1];
  p[0] = 2.0e-8;
  char *filepath;
  filepath = "True_ini";
  printf("Load initial conditions.\n");
  loadInitialPath(xold,filepath);

  char fn1[255];
  sprintf(fn1, "./true_path.dat");
  FILE *fp1 = fopen(fn1,"w");
  printf("Starting integrating.\n");
  for(d=0;d<D;d++){
    fprintf(fp1, "%f\t", xold[d]);
  }
  fprintf(fp1, "\n");
  for(m=0;m<M+1;m++){
    RK4(m,xold,x_mp1,p,10);
    for(d=0;d<D;d++){
      fprintf(fp1, "%f\t", x_mp1[d]);
      xold[d] = x_mp1[d];
    }//end for d
    fprintf(fp1, "\n");
  }//end for m
  fclose(fp1);
  return 0;
}

//forward integrator, 4th order Runge-Kutta
void RK4(int m, double *x_m, double *x_mp1, double *p, int m_step){
  double x1[D],x2[D],x3[D],a1[D],a2[D],a3[D],a4[D];
  int d_temp,step;
  for(d_temp=0;d_temp<D;d_temp++){
    x_mp1[d_temp] = x_m[d_temp];
  }//end for d_temp
  for(step=0;step<m_step;step++){
    for(d_temp=0;d_temp<D;d_temp++){
      a1[d_temp] = dt*f(d_temp,m,x_mp1,p);
      x1[d_temp] = x_mp1[d_temp] + 0.5*a1[d_temp];
    }//end for d_temp
    for(d_temp=0;d_temp<D;d_temp++){
      a2[d_temp] = dt*f(d_temp,m,x1,p);
      x2[d_temp] = x_mp1[d_temp] + 0.5*a2[d_temp];
    }//end for d_temp
    for(d_temp=0;d_temp<D;d_temp++){
      a3[d_temp] = dt*f(d_temp,m,x2,p);
      x3[d_temp] = x_mp1[d_temp] + a3[d_temp];
    }//end for d_temp
    for(d_temp=0;d_temp<D;d_temp++){
      a4[d_temp] = dt*f(d_temp,m,x3,p);
      x_mp1[d_temp] += (a1[d_temp]+2*a2[d_temp]+2*a3[d_temp]+a4[d_temp])/6;
    }//end for d_temp
  }//end for step
}

void loadInitialPath(double *xold, char *fname){
    int d,out;
    double f1;
    FILE *fxold;
    char name[127];
    sprintf(name, "%s.dat", fname);
    fxold = fopen(name,"r");
    if(fxold==NULL){
        printf("File %s for initial path is missing\n",name);
        exit(1);
    }
    printf("Loading initial path from %s\n", name);
    for(d=0;d<D;d++){
          out = fscanf(fxold,"%lf",&f1);
          xold[d] = f1;
    }
    fclose(fxold);
}
