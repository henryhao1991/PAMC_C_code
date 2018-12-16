// Number of variables
#define D 20

// Number of time steps
#define M 200

// Number of parameters
#define NP 1

// Integration step, dt
#define dt 0.025

// Steps between Measurements
#define nstep 1

// Define annealing parameters beta_max and alpha
#define beta_max 75
#define alpha 1.2

// Define Iteration settings
#define N_ini 1000
#define Bsize 500
#define Nblock 10

// Default step size for adjustment at initialization phase
#define delta_0  1.0;//Default step size
#define deltaF_0  1.0;//Default parameter step size

/////////////////////////////////////////////////////////////////////////////////////////
// Put the Differential Equations (dx_d/dt =  ) here
// Variables are y[0], ..., y[D-1],  Parameters are p[0],...,p[NP-1]
/////////////////////////////////////////////////////////////////////////////////////////

double f(int d, int m, double *y, double *p){
    // Lorenz 96
   return  y[(d+D-1)%D] * (y[(d+1)%D] - y[(d+D-2)%D]) - y[d] + p[0];
}
