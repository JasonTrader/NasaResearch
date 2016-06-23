/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and
 * proprietary rights in and to this software and related documentation.
 * Any use, reproduction, disclosure, or distribution of this software
 * and related documentation without an express license agreement from
 * NVIDIA Corporation is strictly prohibited.
 *
 * Please refer to the applicable NVIDIA end user license agreement (EULA)
 * associated with this source code for terms and conditions that govern
 * your use of this NVIDIA software.
 *
 */


#include "book.h"
#include <stdio.h>

#define PI 3.1415926535897932384
#define mu0 4*PI*1e-7

__global__ void init(double *rod_new, double imax, double ldr, double rlength, int totl) {
  int i = blockIdx.x;
  if (i < totl - 1){
    if(i != 0)
      rod_new[i] = (1-(i*i*ldr*ldr/(3*rlength*rlength)))*3*mu0*imax*i*ldr/(4*PI*rlength*rlength);
  }
}

__global__ void copyToOld(double *rod_new, double *rod_old, int totl){
  int i = blockIdx.x;
  if (i < totl - 1)
    rod_old[i] = rod_new[i];
}

int main( void ) {

  FILE *myfile;
  myfile = fopen("results.txt", "w");
  double imax, rlength, eta, tstep, ldr, conv;
  int numseg;

  printf("%s", "What is your I max? ");
  scanf("%lf", &imax);
  printf("%s", "What is the length of your rod? ");
  scanf("%lf", &rlength);
  printf("%s", "What is eta? ");
  scanf("%lf", &eta);
  printf("%s", "How many segments would you like? ");
  scanf("%d", &numseg);
  ldr = rlength/(numseg+1);
  double bound = 0.5*ldr*ldr*mu0/eta;
  printf("%s%lf%s", "What time step would you like? (must be less than ", bound, " ) ");
  scanf("%lf", &tstep);
  printf("%s", "What is the threshold for your convergence? ");
  scanf("%lf", &conv);

  //initialize
  double *rod_new = new double[numseg+2];
  double *rod_old = new double[numseg+2];
  *rod_new = 0;
  *(rod_new + numseg + 1) = mu0*imax/(2*PI*rlength);
  *rod_old = 0;
  *(rod_old + numseg + 1) = mu0*imax/(2*PI*rlength);

  double *dev_old, *dev_new;

  // allocate the memory on the GPU
  HANDLE_ERROR( cudaMalloc( (void**)&dev_old, (numseg+2) * sizeof(double) ) );
  HANDLE_ERROR( cudaMalloc( (void**)&dev_new, (numseg+2) * sizeof(double) ) );

  // fill the arrays 'old' and 'new'
  init<<<numseg+1,1>>>(dev_new, imax, ldr, rlength, numseg+2);

  // copy the array 'c' back from the GPU to the CPU
  HANDLE_ERROR( cudaMemcpy( rod_new, dev_new, (numseg+2) * sizeof(double), cudaMemcpyDeviceToHost ) );

  // output x values
  for (int i=0; i<numseg+2; i++) {
      fprintf( myfile, "%lf ", i*ldr );
  }

  for (int i=0; i<numseg+2; i++) {
      fprintf( myfile, "%lf ", *(rod_new+i) );
  }

  // free the memory allocated on the GPU
  HANDLE_ERROR( cudaFree( dev_old ) );
  HANDLE_ERROR( cudaFree( dev_new ) );

  fclose(myfile);

  return 0;
}
