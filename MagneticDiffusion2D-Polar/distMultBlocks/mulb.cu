#include <stdio.h>
#include <time.h>
using namespace std;

#define PI 3.1415926535897932384
#define mu0 4*PI*1e-7
#define block_i 16
#define block_j 16

//grid will be r driven meaning grid(r,z) = grid[r*zMax + z]

__global__ void init(double *grid, double Il, double dI, double ldr, double rlength, int rseg, int zseg){
  int j = blockIdx.y *block_j + threadIdx.y;
  int i = blockIdx.x *block_i + threadIdx.x;
  if(j<(rseg+2) && i < (zseg+2))
    grid[j*(zseg+2)+i] = (1-(j*j*ldr*ldr/(3*rlength*rlength)))*3*mu0*(Il + i*dI)*j*ldr/(4*PI*rlength*rlength);
}

__global__ void writeOut(double *grid, double *tempGrid, int rseg, int zseg){
  int j = blockIdx.y *block_j + threadIdx.y;
  int i = blockIdx.x *block_i + threadIdx.x;
  if(i<(zseg+2) && j<(rseg+2)){
    grid[j*(zseg+2) + i] = tempGrid[j*(zseg+2) + i];
  }
}

__global__ void run(double *grid, double *tempGrid, double r_aug, double z_aug, int rseg, int zseg){
  int j = blockIdx.y *block_j + threadIdx.y;
  int i = blockIdx.x *block_i + threadIdx.x;
  int sharedPos = threadIdx.y*(block_i+2) + threadIdx.x;

  extern __shared__ double grid_old_s[];

  if(j<(rseg+2) && i<(zseg+2)){
    grid_old_s[sharedPos] = grid[j*(zseg+2) + i];
  }



  __syncthreads();

  if(i<(zseg+1) && i != 0 && j<(rseg+1) && j!=0){
    if(threadIdx.x != 0 && threadIdx.x != block_i + 1 && threadIdx.y != 0 && threadIdx.y != block_j + 1){
      if(j==1){
        tempGrid[j*(zseg+2) + i] += r_aug*(2*grid_old_s[sharedPos+block_i+2] - 4*grid_old_s[sharedPos]) +
        z_aug * (grid_old_s[sharedPos+1] - 2*grid_old_s[sharedPos] + grid_old_s[sharedPos-1]);
      }
      else{
        tempGrid[j*(zseg+2) + i] += r_aug*((1+(1/(2*j)))*grid_old_s[sharedPos+block_i+2] + (-2-(1/(j*j)))*grid_old_s[sharedPos] + (1-(1/(2*j)))*grid_old_s[sharedPos-(block_i+2)])
        +z_aug*(grid_old_s[sharedPos+1] - 2*grid_old_s[sharedPos] + grid_old_s[sharedPos-1]);
      }
    }
  }
/*
  if(i<(zseg+1) && i != 0 && j<(rseg+1) && j!=0){
    if(threadIdx.x != 0 && threadIdx.x != block_i + 1 && threadIdx.y != 0 && threadIdx.y != block_j + 1){
      if(j==1){
        tempGrid[j*(zseg+2) + i] += r_aug*(2*grid[(j+1)*(zseg+2) + i] - 4*grid[j*(zseg+2) + i]) +
        z_aug * (grid[j*(zseg+2) + i + 1] - 2*grid[j*(zseg+2) + i] + grid[j*(zseg+2) + i - 1]);
      }
      else{
        tempGrid[j*(zseg+2) + i] += r_aug*((1+(1/(2*j)))*grid[(j+1)*(zseg+2) + i] + (-2-(1/(j*j)))*grid[j*(zseg+2) + i] + (1-(1/(2*j)))*grid[(j-1)*(zseg+2) + i])
        +z_aug*(grid[j*(zseg+2) + i + 1] - 2*grid[j*(zseg+2) + i] + grid[j*(zseg+2) + i -1]);
      }
    }
  }*/
}

int main(){
  double Il, Ir, rlength, eta, tstep, ldr, ldz, tottime, zlength;
  int rseg, zseg;
  //printf("What is your left I? ");
  scanf("%lf", &Il);
  //printf("What is your right I? ");
  scanf("%lf", &Ir);
  //printf("What is the radius of your rod? ");
  scanf("%lf", &rlength);
  //printf("What is the length of your rod? ");
  scanf("%lf", &zlength);
  //printf("What is eta? ");
  scanf("%lf", &eta);
  //printf("How many segments would you like per radius? ");
  scanf("%d", &rseg);
  //printf("How many segments would you like per length? ");
  scanf("%d", &zseg);
  ldr = rlength/(rseg+1);
  ldz = zlength/(zseg+1);
  double smallest = ldr;
  if(ldz < ldr)
    smallest = ldz;
  tstep = 0.125*smallest*smallest*mu0/eta;
  //printf("How long would you like to run? ");
  scanf("%lf", &tottime);

  double *h_grid, *d_grid, *d_tempGrid;
  size_t grid_size = (rseg + 2)*(zseg+2) * sizeof(double);
  h_grid = (double*)malloc(grid_size);
  cudaMalloc(&d_grid, grid_size);
  cudaMalloc(&d_tempGrid, grid_size);

  double dI = (Ir - Il) / (zseg+2);

  int init_block_i, init_block_j;
  init_block_i = 1+((zseg + 2 - 1)/block_i);
  init_block_j = 1+((rseg+2-1)/block_j);
  dim3 initBlockDim(init_block_i, init_block_j);
  dim3 initThreadDim(block_i, block_j);

  init<<<initBlockDim,initThreadDim>>>(d_grid, Il, dI, ldr, rlength, rseg, zseg);

  cudaMemcpy(h_grid, d_grid, grid_size, cudaMemcpyDeviceToHost);

  FILE *myfile;
  myfile = fopen("init.txt", "w");

  long int i;

  for(i = 0; i< zseg+1; i++)
    fprintf(myfile, "%lf ", i*ldz);
  fprintf(myfile, "%lf\n", i*ldz);

  for(i = 0; i< rseg+1; i++)
    fprintf(myfile, "%lf ", i*ldr);
  fprintf(myfile, "%lf\n", i*ldr);

  for(i = 0; i< (rseg + 2)*(zseg+2); i++){
    if(i%(zseg+2)==zseg+1)
      fprintf(myfile, "%lf\n", h_grid[i]);
    else
      fprintf(myfile, "%lf ", h_grid[i]);

  }

  fclose(myfile);



  double r_aug = eta*tstep/(mu0*ldr*ldr);
  double z_aug = eta*tstep/(mu0*ldz*ldz);

  int run_block_i, run_block_j ;
  run_block_i = 1+((zseg-1)/block_i);
  run_block_j = 1+((rseg-1)/block_j);
  dim3 blockDim(run_block_i, run_block_j);
  //printf("%d\n", run_block_i);
  dim3 threadDim(block_i + 2, block_j + 2);

  cudaMemcpy(d_tempGrid, d_grid, grid_size, cudaMemcpyDeviceToDevice);


  clock_t begin, end;
  double time_spent;
  begin = clock();

  //run
  long int step = 0;
  while((step*tstep) < tottime){
    run<<<blockDim,threadDim, (block_i+2)*(block_j+2)*sizeof(double)>>>(d_grid, d_tempGrid, r_aug, z_aug, rseg, zseg);
    //cudaDeviceSynchronize();
    cudaMemcpy(d_grid, d_tempGrid, grid_size, cudaMemcpyDeviceToDevice);
    //cudaDeviceSynchronize();
    step++;
  }
  cudaDeviceSynchronize();

  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  cudaMemcpy(h_grid, d_grid, grid_size, cudaMemcpyDeviceToHost);
  printf("\nSteps: %ld\n", step);

  myfile = fopen("res.txt", "w");

  for(i = 0; i< zseg+1; i++)
    fprintf(myfile, "%lf ", i*ldz);
  fprintf(myfile, "%lf\n", i*ldz);

  for(i = 0; i< rseg+1; i++)
    fprintf(myfile, "%lf ", i*ldr);
  fprintf(myfile, "%lf\n", i*ldr);

  for(i = 0; i< (rseg + 2)*(zseg+2); i++){
    if(i%(zseg+2)==zseg+1)
      fprintf(myfile, "%lf\n", h_grid[i]);
    else
      fprintf(myfile, "%lf ", h_grid[i]);

  }

  fclose(myfile);


  free(h_grid);
  cudaFree(d_grid);
  cudaFree(d_tempGrid);



  printf("\n------------------------------------\nExecution took: %lf sec\n", time_spent);

  return 0;
}
