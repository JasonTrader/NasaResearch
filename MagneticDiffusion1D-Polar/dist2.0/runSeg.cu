#include <stdio.h>
#include <time.h>

#define PI 3.1415926535897932384
#define mu0 4*PI*1e-7
#define threadsPerBlock 256

__global__ void init(double *rod_new, double imax, double ldr, double rlength, int numseg){
  int i = blockIdx.x*threadsPerBlock + threadIdx.x;
  if(i< numseg+2)
    rod_new[i] = (1-(i*i*ldr*ldr/(3*rlength*rlength)))*3*mu0*imax*i*ldr/(4*PI*rlength*rlength);
}

__global__ void run(double *rod_new, double aug, int numseg){
  int bi = blockIdx.x;
  int ti = threadIdx.x;
  int i = bi*threadsPerBlock + ti;
  extern __shared__ double rod_old_s[];
  if(i< (numseg+2))
    rod_old_s[ti] = rod_new[i];
  __syncthreads();
  /*
  if(bi == 2 && ti == 0){
    printf("%lf\t%d\n", rod_old_s[0],bi);
  }
  if(bi == 1 && ti == 0){
    printf("%lf\t%d\n", rod_old_s[3],bi);
  }
  */
  if(i<(numseg+1) && i != 0 && ti != 0 && ti != (threadsPerBlock + 1)){
    if(i==1)
      rod_new[1]+= aug*(2*rod_old_s[2] - 4*rod_old_s[1]);
    else
      rod_new[i] += aug*((1+(1/(2*i)))*rod_old_s[ti+1] + (-2-(1/(i*i)))*rod_old_s[ti] + (1-(1/(2*i)))*rod_old_s[ti-1]);
  }
}

int main(){
  FILE *myfile;
  myfile = fopen("segResults.txt", "w");
  double imax, rlength, eta, tstep, ldr, tottime;
  int numseg;
  //printf("What is your I max? ");
  scanf("%lf", &imax);
  //printf("What is the length of your rod? ");
  scanf("%lf", &rlength);
  //printf("What is eta? ");
  scanf("%lf", &eta);
  //printf("How many segments would you like? ");
  scanf("%d", &numseg);
  ldr = rlength/(numseg+1);
  tstep = 0.25*ldr*ldr*mu0/eta;
  //printf("How long would you like to run? ");
  scanf("%lf", &tottime);

  double *h_rod, *d_rod;
  size_t rod_size = (numseg + 2) * sizeof(double);
  h_rod = (double*)malloc(rod_size);
  cudaMalloc(&d_rod, rod_size);

//The total number of points passed is numseg + 2, which includes ghost points
  int numInitBlocks = 1 + (((numseg + 2) - 1)/threadsPerBlock);

  init<<<numInitBlocks, threadsPerBlock>>>(d_rod, imax, ldr, rlength, numseg);

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
  //long int total_steps = tottime / tstep;

//Every block will get its own ghost cell value, so don't need the +2 for numseg
  int numBlocks = 1 + ((numseg - 1)/threadsPerBlock);
  clock_t begin, end;
  double time_spent;
  begin = clock();

  //run
  long int steps = 0;
  while(steps*tstep< tottime){
    run<<<numBlocks, threadsPerBlock + 2, (threadsPerBlock+2)*sizeof(double)>>>(d_rod, aug, numseg);
    steps++;
  }
  cudaDeviceSynchronize();

  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  printf("\nSteps: %ld\n", steps);
  cudaMemcpy(h_rod, d_rod, rod_size, cudaMemcpyDeviceToHost);

  for(out = 0; out<numseg+1; out++){
    fprintf( myfile, "%lf ", *(h_rod+out) );
  }
  fprintf( myfile, "%lf\n", *(h_rod+out) );

  fprintf(myfile, "STOP\n");
  fclose(myfile);

  cudaFree(d_rod);
  free(h_rod);



  printf("\n------------------------------------\nExecution took: %lf sec\n", time_spent);

  return 0;
}
