#ifndef FUNC_H
#define FUNC_H

#define PI 3.1415926535897932384
#define mu0 4*PI*1e-7

__global__ void init(double *rod_new, double imax, double ldr, double rlength, int totl) {
  int i = blockIdx.x;
  if (i < totl)
    rod_new[i] = (1-(i*i*ldr*ldr/(3*rlength*rlength)))*3*mu0*imax*i*ldr/(4*PI*rlength*rlength);
}

__global__ void copyToOld(double *rod_new, double *rod_old, int totl){
  int i = blockIdx.x;
  if (i < totl)
    rod_old[i] = rod_new[i];
}

__global__ void update(double *rod_new, double *rod_old, int totl, bool *conv, double thresh, double aug){
  int i = blockIdx.x + 1;
  if (i < totl){
    if(i == 1)
      rod_new[1]+= aug*(2*rod_old[2] - 4*rod_old[1]);
    else
      rod_new[i] += aug*((1+(1/(2*i)))*rod_old[i+1] + (-2-(1/(i*i)))*rod_old[i] + (1-(1/(2*i)))*rod_old[i-1]);

    double diff = rod_new[i] - rod_old[i];
    if (diff < 0)
      diff *= -1;
    conv[i-1] = (diff < thresh);
  }
}

bool converge(bool *conv, int length){
  for(int i = 0; i< length; i++){
    if(!conv[i])
      return false;
  }
  return true;
}

#endif
