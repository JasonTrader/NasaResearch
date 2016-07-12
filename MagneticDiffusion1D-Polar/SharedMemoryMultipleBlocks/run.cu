#include <iostream>
#include <stdio.h>
#include <time.h>
#include <iomanip>
using namespace std;

#define PI 3.1415926535897932384
#define mu0 4*PI*1e-7
#define threadsPerBlock 1000

__global__ void setGhostPoints(double *rod_new, double *ghost, int numBlocks, int numseg, long int steps){
  int sec;
  threadIdx.x == numBlocks - 1 ? sec = numseg + 1 : sec = (threadIdx.x+1) * threadsPerBlock + 1;
  ghost[2*threadIdx.x] = rod_new[threadIdx.x*threadsPerBlock];
  ghost[2*threadIdx.x+1] = rod_new[sec];
}

__global__ void init(double *rod_new, double imax, double ldr, double rlength, int numseg, int numBlocks){
  int i = blockIdx.x*threadsPerBlock + threadIdx.x + 1;
  rod_new[i] = (1-(i*i*ldr*ldr/(3*rlength*rlength)))*3*mu0*imax*i*ldr/(4*PI*rlength*rlength);
}

__global__ void initBound(double *rod_new, double imax, double ldr, double rlength, int numseg){
  rod_new[numseg+1] = mu0*imax/(2*PI*rlength);
  rod_new[0] = 0;
}

__global__ void run(double *ghost, double *rod_new, double aug, int numseg, int numBlocks, long int steps){
  int bi = blockIdx.x;
  int ti = threadIdx.x;
  int i = bi*threadsPerBlock + ti + 1;
  int threadsNeeded;
  bi == (numBlocks - 1) ? threadsNeeded = (numseg - (bi*threadsPerBlock)) : threadsNeeded =  threadsPerBlock;
  extern __shared__ double rod_old_s[];
  double ghost_left = ghost[2*bi];
  double ghost_right = ghost[2*bi+1];
  rod_old_s[ti] = rod_new[i];
  __syncthreads();

  if(threadsNeeded == 1 && ti == 0)
    rod_new[i] += aug*((1+(1/(2*i)))*ghost_right + (-2-(1/(i*i)))*rod_old_s[ti] + (1-(1/(2*i)))*ghost_left);
  else if(i==1)
    rod_new[1]+= aug*(2*rod_old_s[ti+1] - 4*rod_old_s[ti]);
  else if(ti == 0)
    rod_new[i] += aug*((1+(1/(2*i)))*rod_old_s[ti+1] + (-2-(1/(i*i)))*rod_old_s[ti] + (1-(1/(2*i)))*ghost_left);
  else if(ti == threadsNeeded - 1)
    rod_new[i] += aug*((1+(1/(2*i)))*ghost_right + (-2-(1/(i*i)))*rod_old_s[ti] + (1-(1/(2*i)))*rod_old_s[ti-1]);
  else if(i<(numseg + 1))
    rod_new[i] += aug*((1+(1/(2*i)))*rod_old_s[ti+1] + (-2-(1/(i*i)))*rod_old_s[ti] + (1-(1/(2*i)))*rod_old_s[ti-1]);
}

int main(){
  FILE *myfile;
  myfile = fopen("segResults.txt", "w");
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

  double *h_rod, *d_rod, *d_ghost;
  size_t rod_size = (numseg + 2) * sizeof(double);
  h_rod = (double*)malloc(rod_size);
  cudaMalloc(&d_rod, rod_size);

  int numBlocks = (numseg + threadsPerBlock -1)/threadsPerBlock;
  size_t ghost_size = numBlocks * 2 * sizeof(double);
  cudaMalloc(&d_ghost, ghost_size);

  initBound<<<1,1>>>(d_rod, imax, ldr, rlength, numseg);
  init<<<numBlocks,threadsPerBlock>>>(d_rod, imax, ldr, rlength, numseg, numBlocks);

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
  printf("\nSteps: %ld\tTstep: %le\n", total_steps, tstep);


  clock_t begin, end;
  double time_spent;
  begin = clock();

    //run
  long int steps = 0;
  while(steps< total_steps){
    setGhostPoints<<<1,numBlocks>>>(d_rod, d_ghost, numBlocks, numseg, steps);
    run<<<numBlocks, threadsPerBlock, threadsPerBlock*sizeof(double)>>>(d_ghost, d_rod, aug, numseg, numBlocks, steps);
    steps++;
  }

  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  cudaMemcpy(h_rod, d_rod, rod_size, cudaMemcpyDeviceToHost);

  for(out = 0; out<numseg+1; out++){
    fprintf( myfile, "%lf ", *(h_rod+out) );
  }
  fprintf( myfile, "%lf\n", *(h_rod+out) );

  fprintf(myfile, "STOP\n");
  fclose(myfile);

  cudaFree(d_rod);
  cudaFree(d_ghost);
  free(h_rod);



  cout << "\n------------------------------------\nExecution took: "<<  time_spent << " sec\n";

  return 0;
}
