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

__global__ void seeOld(double *rod_old, int totl){
  int i = blockIdx.x;
  if (i < totl)
    printf("\n%d: %lf", i, rod_old[i]);
}

__device__ double get_next(double *rod_old, int totl, int num_steps, double aug, int pos, double *temp){
  bool first_complete = false;
  bool set = false;
  for(int step = (num_steps-1); step>-1; step--){
    for(int i = (pos-step); i<(pos+step+1); i++){
      int share_mem_pos = i-pos+num_steps-1+((pos-1)*((2*num_steps)-1));
      if(!first_complete){
        set = true;
        if(i>1 && i<totl){
          temp[share_mem_pos] = rod_old[i] + aug*((1+(1/(2*i)))*rod_old[i+1] + (-2-(1/(i*i)))*rod_old[i] + (1-(1/(2*i)))*rod_old[i-1]);
        }
        else if(i==1){
          temp[share_mem_pos] = rod_old[i] + aug*(2*rod_old[2] - 4*rod_old[1]);
        }
      }
      else{
        if(i>1 && i<totl-1){
          temp[share_mem_pos] += aug*((1+(1/(2*i)))*temp[share_mem_pos+1] + (-2-(1/(i*i)))*temp[share_mem_pos] + (1-(1/(2*i)))*temp[share_mem_pos-1]);
        }
        else if(i==1){
          temp[share_mem_pos] += aug*(2*temp[share_mem_pos+1] - 4*temp[share_mem_pos]);
        }
        else if(i==(totl-1)){
          temp[share_mem_pos] += aug*((1+(1/(2*i)))*rod_old[i+1] + (-2-(1/(i*i)))*temp[share_mem_pos] + (1-(1/(2*i)))*temp[share_mem_pos-1]);
        }
      }
    }
    first_complete = set;
  }
  //printf("\n%lf", temp[(pos-1)*num_steps]);
  return temp[(pos-1)*((2*num_steps)-1)+num_steps-1];
}

__global__ void update(double *rod_new, double *rod_old, int totl, bool *conv, double thresh, double aug, int num_steps, double *temp){
  int i = blockIdx.x + 1;
  if (i < totl){
    rod_new[i] = get_next(rod_old, totl, num_steps, aug, i, temp);
    double diff = rod_new[i] - rod_old[i];
    if (diff < 0)
      diff *= -1;
    conv[i-1] = (diff < (thresh*num_steps));
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
