#include <stdio.h>
#include <time.h>
using namespace std;

#define PI 3.1415926535897932384
#define mu0 4*PI*1e-7
#define threadsPerBlock 1024

//grid will be r driven meaning grid(r,z) = grid[r*zMax + z]

__global__ void init(double *grid, double Il, double dI, double ldr, double rlength, int grid_size, int zseg){
  int rem = grid_size%threadsPerBlock;
  int divi = grid_size/threadsPerBlock;
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
    long int rrow = i/(zseg+2);
    long int zcol = i%(zseg+2);
    grid[i] = (1-(rrow*rrow*ldr*ldr/(3*rlength*rlength)))*3*mu0*(Il + zcol*dI)*rrow*ldr/(4*PI*rlength*rlength);
  }
}

__global__ void run(double *rod_new, double r_aug, double z_aug, long int maxSteps, int grid_size, int rseg, int zseg){
  int rem = grid_size%threadsPerBlock;
  int divi = grid_size/threadsPerBlock;
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
  extern __shared__ double grid_new_s[];
  extern __shared__ double grid_old_s[];

  for(int i = start; i<fin; i++){
    grid_new_s[i] = rod_new[i];
  }
  __syncthreads();

  while(steps<maxSteps){
    for(int i = start; i<fin; i++){
      grid_old_s[i] = grid_new_s[i];
    }
    __syncthreads();
    for(int i = start; i<fin; i++){
      int r = i/(zseg+2);
      int z = i%(zseg+2);

      if(z != 0 && z != zseg+1){
        if(r!= 0 && r!= rseg+1){
          if(r==1){
            grid_new_s[i] += r_aug*(2*grid_old_s[i+(zseg+2)] - 4*grid_old_s[i]) +
              z_aug * (grid_old_s[i+1] - 2*grid_old_s[i] + grid_old_s[i-1]);
          }
          else{
            grid_new_s[i] += r_aug*((1+(1/(2*i)))*grid_old_s[i+(zseg+2)] + (-2-(1/(i*i)))*grid_old_s[i] + (1-(1/(2*i)))*grid_old_s[i-(zseg+2)])
              +z_aug*(grid_old_s[i+1] - 2*grid_old_s[i] + grid_old_s[i-1]);
          }
        }
      }
    }
    steps++;
    __syncthreads();
  }
  for(int i = start; i<fin; i++){
    rod_new[i] = grid_new_s[i];
  }
}

int main(){
  double Il, Ir, rlength, eta, tstep, ldr, ldz, tottime, zlength;
  int rseg, zseg;
  printf("What is your left I? ");
  scanf("%lf", &Il);
  printf("What is your right I? ");
  scanf("%lf", &Ir);
  printf("What is the radius of your rod? ");
  scanf("%lf", &rlength);
  printf("What is the length of your rod? ");
  scanf("%lf", &zlength);
  printf("What is eta? ");
  scanf("%lf", &eta);
  printf("How many segments would you like per radius? ");
  scanf("%d", &rseg);
  printf("How many segments would you like per length? ");
  scanf("%d", &zseg);
  ldr = rlength/(rseg+1);
  ldz = zlength/(zseg+1);
  double smallest = ldr;
  if(ldz < ldr)
    smallest = ldz;
  tstep = 0.125*smallest*smallest*mu0/eta;
  printf("How long would you like to run? ");
  scanf("%lf", &tottime);

  double *h_grid, *d_grid;
  size_t grid_size = (rseg + 2)*(zseg+2) * sizeof(double);
  h_grid = (double*)malloc(grid_size);
  cudaMalloc(&d_grid, grid_size);

  double dI = (Ir - Il) / (zseg+2);

  init<<<1,threadsPerBlock>>>(d_grid, Il, dI, ldr, rlength, (rseg + 2)*(zseg+2), zseg);

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
  long int total_steps = tottime / tstep;
  printf("\nSteps: %ld\n", total_steps);


  clock_t begin, end;
  double time_spent;
  begin = clock();

  //run
  printf("Called run\n");
  run<<<1,threadsPerBlock, (rseg + 2)*(zseg+2)*sizeof(double)>>>(d_grid, r_aug, z_aug, total_steps, (rseg + 2)*(zseg+2), rseg, zseg);
  cudaDeviceSynchronize();

  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  cudaMemcpy(h_grid, d_grid, grid_size, cudaMemcpyDeviceToHost);

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



  printf("\n------------------------------------\nExecution took: %lf sec\n", time_spent);

  return 0;
}
