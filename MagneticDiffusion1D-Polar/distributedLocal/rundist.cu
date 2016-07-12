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
#include "func.h"
#include <stdio.h>
#include <time.h>

#define PI 3.1415926535897932384
#define mu0 4*PI*1e-7

int main( void ) {



  FILE *myfile;
  myfile = fopen("results.txt", "w");
  double imax, rlength, eta, tstep, ldr, tottime;
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
  printf("%s", "How long would you like to run?  \n");
  scanf("%lf", &tottime);

  //initialize
  double *rod_new = new double[numseg+2];

  double *rold, *rnew;
  bool *ready;
  int *timestep;

  // allocate the memory on the GPU
  HANDLE_ERROR( cudaMalloc( (void**)&rold, (numseg+2) * sizeof(double) ) );
  HANDLE_ERROR( cudaMalloc( (void**)&rnew, (numseg+2) * sizeof(double) ) );
  HANDLE_ERROR( cudaMalloc( (void**)&ready, (numseg+2) * sizeof(bool) ) );
  HANDLE_ERROR( cudaMalloc( (void**)&timestep, (numseg+2) * sizeof(int) ) );

  // fill the array 'new'
  init<<<numseg+2,1>>>(rnew, rold, timestep, ready, imax, ldr, rlength, numseg+2);

  // copy data on device from 'rnew' to 'rod_new'
  HANDLE_ERROR( cudaMemcpy( rod_new, rnew, (numseg+2) * sizeof(double), cudaMemcpyDeviceToHost ) );


  int out;
  // output r values
  for (out = 0; out<numseg+1; out++) {
      fprintf( myfile, "%lf ", out*ldr );
  }
  fprintf( myfile, "%lf\n", out*ldr );

  double aug = eta*tstep/(mu0*ldr*ldr);

  //output initial conditions
  for (out=0; out<numseg+1; out++) {
      fprintf( myfile, "%lf ", *(rod_new+out) );
  }
  fprintf( myfile, "%lf\n", *(rod_new+out) );

  int steps = tottime / tstep;

  clock_t begin, end;
  double time_spent;
  begin = clock();
  //run
  run<<<numseg,1>>>(rnew, rold, numseg+1, aug, steps, ready, timestep);
  HANDLE_ERROR( cudaMemcpy( rod_new, rnew, (numseg+2) * sizeof(double), cudaMemcpyDeviceToHost ) );

  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  // free the memory allocated on the GPU
  HANDLE_ERROR( cudaFree( rold ) );
  HANDLE_ERROR( cudaFree( rnew ) );
  HANDLE_ERROR( cudaFree( ready ) );
  HANDLE_ERROR( cudaFree( timestep ) );

//output final values
  for (out=0; out<numseg+1; out++) {
      fprintf( myfile, "%lf ", *(rod_new+out) );
  }
  fprintf( myfile, "%lf\n", *(rod_new+out) );

  fprintf(myfile, "STOP\n");
  fclose(myfile);



  printf("\n------------------------------------\n");
  printf("Execution took: %lf sec\n", time_spent);
  return 0;
}
