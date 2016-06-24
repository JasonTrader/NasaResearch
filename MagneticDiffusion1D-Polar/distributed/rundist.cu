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

  clock_t begin, end;
  double time_spent;
  begin = clock();

  FILE *myfile;
  myfile = fopen("results.txt", "w");
  double imax, rlength, eta, tstep, ldr, thresh;
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
  scanf("%lf", &thresh);

  //initialize
  double *rod_new = new double[numseg+2];
  bool *conv = new bool[numseg];

  double *dev_old, *dev_new;
  bool *dev_conv;

  // allocate the memory on the GPU
  HANDLE_ERROR( cudaMalloc( (void**)&dev_old, (numseg+2) * sizeof(double) ) );
  HANDLE_ERROR( cudaMalloc( (void**)&dev_new, (numseg+2) * sizeof(double) ) );
  HANDLE_ERROR( cudaMalloc( (void**)&dev_conv, numseg * sizeof(bool) ) );

  // fill the array 'new'
  init<<<numseg+2,1>>>(dev_new, imax, ldr, rlength, numseg+2);

  // copy data on device from 'dev_new' to 'rod_new'
  HANDLE_ERROR( cudaMemcpy( rod_new, dev_new, (numseg+2) * sizeof(double), cudaMemcpyDeviceToHost ) );


  int out;
  // output r values
  for (out = 0; out<numseg+1; out++) {
      fprintf( myfile, "%lf ", out*ldr );
  }
  fprintf( myfile, "%lf\n", out*ldr );

  double aug = eta*tstep/(mu0*ldr*ldr);
  int tcount = 0;

  do{
    printf("\ntimestep %d", tcount);
    if(tcount%10000==0){
      HANDLE_ERROR( cudaMemcpy( rod_new + 1, dev_new + 1, numseg * sizeof(double), cudaMemcpyDeviceToHost ) );
      for (out=0; out<numseg+1; out++) {
          fprintf( myfile, "%lf ", *(rod_new+out) );
      }
      fprintf( myfile, "%lf\n", *(rod_new+out) );
    }
    tcount++;

    //copy new to old
    copyToOld<<<numseg+2,1>>>(dev_new, dev_old, numseg+2);



    //update
    update<<<numseg,1>>>(dev_new, dev_old, numseg+1, dev_conv, thresh, aug);
    HANDLE_ERROR( cudaMemcpy( conv, dev_conv, numseg, cudaMemcpyDeviceToHost ) );

  } while(!converge(conv, numseg));

  // free the memory allocated on the GPU
  HANDLE_ERROR( cudaFree( dev_old ) );
  HANDLE_ERROR( cudaFree( dev_new ) );
  HANDLE_ERROR( cudaFree( dev_conv ) );

  fprintf(myfile, "STOP\n");
  fclose(myfile);

  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  printf("\n------------------------------------\n");
  printf("Execution took: %lf sec\n", time_spent);
  return 0;
}
