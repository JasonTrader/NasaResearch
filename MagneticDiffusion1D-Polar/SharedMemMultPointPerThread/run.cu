#include <iostream>
#include <stdio.h>
#include <time.h>
using namespace std;

#define PI 3.1415926535897932384
#define mu0 4*PI*1e-7
#define threadsPerBlock 1024

__global__ void init(double *rod_new, double imax, double ldr, double rlength, int rod_size){
  int rem = rod_size%threadsPerBlock;
  int divi = rod_size/threadsPerBlock;
  int start, fin;
  if(threadIdx.x<rem){
    start = threadIdx.x*(divi+1);
    fin = start + divi + 1;
  }
  else{
    start = threadIdx.x*divi + rem;
    fin = start + divi;
  }
  for(int i = start; i<fin; i++){
    rod_new[i] = (1-(i*i*ldr*ldr/(3*rlength*rlength)))*3*mu0*imax*i*ldr/(4*PI*rlength*rlength);
  }
}

__global__ void run(double *rod_new, double aug, long int maxSteps, int rod_size){
  int rem = rod_size%threadsPerBlock;
  int divi = rod_size/threadsPerBlock;
  int start, fin;
  if(threadIdx.x<rem){
    start = threadIdx.x*(divi+1);
    fin = start + divi + 1;
  }
  else{
    start = threadIdx.x*divi + rem;
    fin = start + divi;
  }
  long int steps = 0;
  extern __shared__ double rod_new_s[];
  extern __shared__ double rod_old_s[];
  for(int i = start; i<fin; i++){
    rod_new_s[i] = rod_new[i];
  }
  __syncthreads();

  while(steps<maxSteps){
    for(int i = start; i<fin; i++){
      rod_old_s[i] = rod_new_s[i];
    }
    __syncthreads();
    for(int i = start; i<fin; i++){
      if(i==1)
        rod_new_s[1]+= aug*(2*rod_old_s[2] - 4*rod_old_s[1]);
      else if(i<(rod_size - 1) && i > 1)
        rod_new_s[i] += aug*((1+(1/(2*i)))*rod_old_s[i+1] + (-2-(1/(i*i)))*rod_old_s[i] + (1-(1/(2*i)))*rod_old_s[i-1]);
      }
    steps++;
    __syncthreads();
  }
  for(int i = start; i<fin; i++){
    rod_new[i] = rod_new_s[i];
  }
}

int main(){
  FILE *myfile;
  myfile = fopen("results.txt", "w");
  double imax, rlength, eta, tstep, ldr, tottime;
  int numseg;
  printf("What is your I max? ");
  scanf("%lf", &imax);
  printf("What is the length of your rod? ");
  scanf("%lf", &rlength);
  printf("What is eta? ");
  scanf("%lf", &eta);
  printf("How many segments would you like? ");
  scanf("%d", &numseg);
  ldr = rlength/(numseg+1);
  tstep = 0.25*ldr*ldr*mu0/eta;
  printf("How long would you like to run? ");
  scanf("%lf", &tottime);

  double *h_rod, *d_rod;
  size_t rod_size = (numseg + 2) * sizeof(double);
  h_rod = (double*)malloc(rod_size);
  cudaMalloc(&d_rod, rod_size);

  init<<<1,threadsPerBlock>>>(d_rod, imax, ldr, rlength, numseg + 2);

  int out;
  //output r values
  for(out = 0; out<numseg+1; out++){
    fprintf( myfile, "%lf ", out*ldr );
  }
  fprintf( myfile, "%lf\n", out*ldr );

  cudaMemcpy(h_rod, d_rod, rod_size, cudaMemcpyDeviceToHost);

  for(out = 0; out<numseg+1; out++){
    fprintf( myfile, "%lf ", *(h_rod+out) );
  }
  fprintf( myfile, "%lf\n", *(h_rod+out) );

  double aug = eta*tstep/(mu0*ldr*ldr);
  long int total_steps = tottime / tstep;
  printf("\nSteps: %ld\n", total_steps);


  clock_t begin, end;
  double time_spent;
  begin = clock();

  //run
  run<<<1,threadsPerBlock, (numseg+2)*sizeof(double)>>>(d_rod, aug, total_steps, numseg+2);
  cudaDeviceSynchronize();

  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  cudaMemcpy(h_rod, d_rod, rod_size, cudaMemcpyDeviceToHost);

  for(out = 0; out<numseg+1; out++){
    fprintf( myfile, "%lf ", *(h_rod+out) );
  }
  fprintf( myfile, "%lf\n", *(h_rod+out) );

  fprintf(myfile, "STOP\n");
  fclose(myfile);

  free(h_rod);
  cudaFree(d_rod);



  cout << "\n------------------------------------\nExecution took: "<<  time_spent << " sec\n";

  return 0;
}
