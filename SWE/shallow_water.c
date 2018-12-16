#include <math.h>
#include <float.h>

#include "model_swe.h"

////////////////////////////////////////////////////////////////////////
//return the neighor index of x(dx,dy), direction defines as follows:
//                             1
//12 13 14 15                  |
//8  9  10 11   +y         2-------0
//4  5  6  7     |             |
//0  1  2  3     ----+x        3
////////////////////////////////////////////////////////////////////////
double n(int dx, int dy, int dir){
  if(dir==0){
    return dy*xdim+((dx+1)%xdim);
  }
  else if(dir==1){
    return ((dy+1)%ydim)*xdim+dx;
  }
  else if(dir==2){
    return dy*xdim+((dx+xdim-1)%xdim);
  }
  else if(dir==3){
    return ((dy+ydim-1)%ydim)*xdim+dx;
  }
  else{
    printf("Direction can only be 0,1,2,3. exit.");
    exit(1);
  }

double eta(int dx, int dy, double *x){
  int d = dy*xdim+dx;
  double vdx, udy, pAve_xy;
  vdx = (x[n(dx,dy,0)+DB]-x[d+DB])/lenx;
  udy = (x[n(dx,dy,1)]-x[d])/leny;
  pAve_xy = (x[d+2*DB]+x[n(dx,dy,0)+2*DB]+x[n(dx,dy,1)+2*DB]+x[n((dx+1%xdim),dy,1)+2*DB])/4.0;
  reutrn (vdx-udy)/pAve_xy;
}

double U(int dx, int dy, double *x){
  int d = dy*xdim+dx;
  return x[d]*(x[d+2*DB]+x[n(dx,dy,0)+2*DB])/2.0;
}

double V(int dx, int dy, double *x){
  int d = dy*xdim+dx;
  return x[d+DB]*(x[d+2*DB]+x[n(dx,dy,1)+2*DB])/2.0;
}

double H(int dx, int dy, double *x){
  int d = dy*xdim+dx;
  return x[d+2*DB]+(x[d]*x[d]+x[n(dx,dy,2)]*x[n(dx,dy,2)]+x[d+DB]*x[d+DB]+
    x[n(dx,dy,3)+DB]*x[n(dx,dy,3)+DB])/4.0;
}

double f(int d, int m, double *x, double *p){
  // Shallow water equations
  int dx,dy,field;
  //relevant field, in order of x,0,1,2,3 as in neighor function
  double eta[5];
  double xdt;

  field = (d-(d%(xdim*ydim)))/(xdim*ydim);
  dx = (d%(xdim*ydim))%xdim;
  dy = ((d%(xdim*ydim))-dx)/dx;
  if(field==0){
    xdt = (eta(dx,(dy+ydim-1)%ydim,x)+eta(dx,dy,x))*(V(dx,dy,x)+
      V((dx+1)%xdim,dy,x)+V(dx,(dy+ydim-1)%ydim,x)+V((dx+1)%xdim,(dy+ydim-1)%ydim,x))
      /8.0-(H((dx+1)%xdim,dy,x)-H(dx,dy,x))/lenx;
    xdt += (f0 + be*dy*leny)*(x[dy*xdim+dx+DB]+x[n(dx,dy,0)+DB]+x[n(dx,dy,3)+DB]
      +x[((dy+ydim-1)%ydim)*xdim+((dx+xdim-1)%xdim)+DB]);
    xdt += dissipation*((x[n(dx,dy,0)]+x[n(dx,dy,2)]-2.0*x[d])/lenx/lenx+
      (x[n(dx,dy,1)]+x[n(dx,dy,3)]-2.0*x[d])/leny/leny);
    xdt += -p[0]*x[d];
    xdt += wind*x[d]*cos(2.0*dy*M_PI/ydim);
  }
  else if(field==1){
    xdt = -(eta((dx+xdim-1)%xdim,dy,x)+eta(dx,dy,x))*(U(dx,dy,x)+
      U((dx+xdim-1)%xdim,dy,x)+U(dx,(dy+1)%ydim,x)+U((dx+xdim-1)%xdim,(dy+1)%ydim,x))
      /8.0-(H(dx,(dy+1)%ydim,x)-H(dx,dy,x))/lenx;
      xdt -= (f0 + be*dy*leny)*(x[dy*xdim+dx]+x[n(dx,dy,1)]+x[n(dx,dy,2)]
        +x[((dy+1)%ydim)*xdim+((dx+xdim-1)%xdim)]);
      xdt += dissipation*((x[n(dx,dy,0)+DB]+x[n(dx,dy,2)+DB]-2.0*x[d])/lenx/lenx+
        (x[n(dx,dy,1)+DB]+x[n(dx,dy,3)+DB]-2.0*x[d])/leny/leny);
      xdt += -p[0]*x[d];
  }
  else if(field==2){
    xdt = -(U(dx,dy,x)-U((dx+xdim-1)%xdim,dy,x))/lenx
      -(V(dx,dy,x)-V(dx,(dy+ydim-1)%ydim,x))/leny;
  }
  else{
    printf("Something wrong with index. (Field can only be 0,1,2.)");
    exit(1);
  }
  return xdt;
}
