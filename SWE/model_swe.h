// Number of time steps
#define M 200

// Number of parameters
#define NP 1

// Integration step, dt
#define dt 36.0

// Steps between Measurements
#define nstep 10

// x and y dimensions
#define xdim 16
#define ydim 16

// Number of variables
#define D xdim*ydim*3

//Shallow water specific constants
#define lenx 5.0e4
#define leny 5.0e4
#define DB xdim*ydim
#define f0 5.0e-5
#define be 2.0e-11
#define dissipation 1.0e-4
#define wind 0.00002
#define Amp 1.0e6
#define depth 5100.0

// Define annealing parameters beta_max and alpha
#define beta_max 75
#define alpha 1.2

// Define Iteration settings
#define N_ini 1000
#define Bsize 500
#define Nblock 10

// Default step size for adjustment at initialization phase
#define delta_0 = 1.0//Default step size
#define deltaF_0 = 1.0e-8//Default parameter step size

/////////////////////////////////////////////////////////////////////////////////////////
// Put the Differential Equations (dx_d/dt =  ) here
// Variables are y[0], ..., y[D-1],  Parameters are p[0],...,p[NP-1]
/////////////////////////////////////////////////////////////////////////////////////////

int n(int dx, int dy, int dir);
double eta(int dx, int dy, double *x);
double U(int dx, int dy, double *x);
double V(int dx, int dy, double *x);
double H(int dx, int dy, double *x);
double f(int d, int m, double *x, double *p);
