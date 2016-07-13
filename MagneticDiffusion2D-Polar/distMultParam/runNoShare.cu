#include <stdio.h>
#include <time.h>

#define PI 3.1415926535897932384
#define mu0 4*PI*1e-7
//Threads per block is capped at 1024 for hardware reasons
//In some cases using a smaller number of threads per block will be more efficient
#define threadsPerBlock 1024
//Max grid points is to defined in order to allocate shared memory for each block
#define MaxGridPoints 6144


/*
NOTE:
grid is technically a 3-dimensional variable length array. In order to copy
this 3d array, it needed to be compressed into a single dimension.
Rather than decompressing this and wasting valuable resources all computation is done in line.
i.e. all blocks will write to adjacent memory
the position in memory is determined by three factors: the block numnber, the r
segment value and the z segment value
Each block is given an initial offset. This offset gives the block a unique space inside
of the overall shared array in which to write its data. Each block's grid is
r driven meaning grid(r,z) = grid[r*zMax + z], thus the position of each point in the
overall array is given by grids[offset + r*zMax + z]
*/


//This function will be called to initialize all data in the grid, Each thread will
//calculate the initial condion for a variable number of points and each block will have a different
//set of initial conditions all passed to device memory already
__global__ void init(double *grids, double *Ils, double *dIs, double *ldrs, double *rlengths, long int *grid_sizes, int *zsegs, long int *offsets){
  // get needed data from memory, all memory is block adjacent
  //block index determines which problem you are working on

  //Current on the left side of the grid
  double Il = Ils[blockIdx.x];
  //amount the current increments per segment in the z direction
  //This model assumes linear change
  double dI = dIs[blockIdx.x];
  //Physical distance change between r segments
  double ldr = ldrs[blockIdx.x];
  //total r length
  double rlength = rlengths[blockIdx.x];
  //Total number of grid points
  long int grid_size = grid_sizes[blockIdx.x];
  //number of z segmnents
  int zseg = zsegs[blockIdx.x];
  //Offset of the current grid
  //Because data needs to be copied to the gpu, it is compressed into one dimension
  //thus all grid data lives in the same array and the offset for that grid determines the satrt of this grid
  //This is explained in more detail above
  long int gridStart = offsets[blockIdx.x];

  //Find how many points per thread will be executing

  //remainder tells you how many threads will have one extra point to cover
  int rem = grid_size%threadsPerBlock;
  //division tells you how many points ALL theads will cover
  int divi = grid_size/threadsPerBlock;
  long int start, fin;
  //If threadIdx.x is less than rem the current thread will recieve one extra point to calculate
  if(threadIdx.x<rem){
    start = threadIdx.x*(divi+1);
    fin = start + divi + 1;
  }
  else{
    //rem is added here because at this point a "rem" number of threads will have one extra point to compute
    //overall this is "rem" extra points and this is accounted for by adding rem
    start = threadIdx.x*divi + rem;
    fin = start + divi;
  }

  //get initial conditions for all points per thread

  //Loops over each point a thread takes care of
  //NOTE: If a thread is not supposed to execute anything it will have the same
  //value for "start" and "fin" and will not enter the loop
  for(int i = start; i<fin; i++){
    //Gets current r and z segmenty values
    //Because grid is "r-driven", dividing by zMax = the r value and moding by
    //zMax gves the current z value
    //For more info check out the "note" at the top
    long int r_seg_val = i/(zseg+2);
    long int z_seg_val = i%(zseg+2);
    //Sets initial condition at the given point based on an equation given by Dr. Sankaran
    grids[i + gridStart] = (1-(r_seg_val*r_seg_val*ldr*ldr/(3*rlength*rlength)))*3*mu0*(Il + z_seg_val*dI)*r_seg_val*ldr/(4*PI*rlength*rlength);
  }
}

//This function funs until the number of timesteps needed is completed
//again this function passes data by memory copy in the form of arrays
//each block excercises a different set of input params
__global__ void run(double *grids, double *r_augs, double *z_augs, long int *allMaxSteps, long int *grid_sizes, int *rsegs, int *zsegs, long int *offsets, double *grids_old){
  //get data from initialization arrays
  //block index determines the problem that is being worked on

  //r_aug and z_aug are constants that multiple values need to update to 1 timestep ahead
  //r_aug is used for changes in r and z_aug is used for changes in z
  double r_aug = r_augs[blockIdx.x];
  double z_aug = z_augs[blockIdx.x];
  //max steps is the number of steps needed to be completed for the simulation
  //to be finished
  long int maxSteps = allMaxSteps[blockIdx.x];
  //total grid size for this block's problem
  long int grid_size = grid_sizes[blockIdx.x];
  //number of segments in r and z directions
  int rseg = rsegs[blockIdx.x];
  int zseg = zsegs[blockIdx.x];
  //offset for block's memory locations
  //More info at top
  long int gridStart = offsets[blockIdx.x];

  //find how many points each thread will execute for

  //rem represents the number of threads that will have one more than divi points
  //to take care of
  int rem = grid_size%threadsPerBlock;
  //divi is the amount of points each thread will have to take care of at a minimum
  int divi = grid_size/threadsPerBlock;
  int start, fin;
  if(threadIdx.x<rem){
    start = threadIdx.x*(divi+1);
    fin = start + divi + 1;
  }
  else{
    //rem is added here because at this point a "rem" number of threads will have one extra point to compute
    //overall this is "rem" extra points and this is accounted for by adding rem
    start = threadIdx.x*divi + rem;
    fin = start + divi;
  }
  long int steps = 0;


  while(steps<maxSteps){
    //copy all points from new to old
    for(int i = start; i<fin; i++){
      grids_old[i + gridStart] = grids[i + gridStart];
    }
    //wait for all threads
    __syncthreads();
    for(int i = start; i<fin; i++){
      //because the grid is row driven
      //i / zMax = rvalue
      //i % zMax = zvalue
      int r = i/(zseg+2);
      int z = i%(zseg+2);

      //leave boundary conditions alone as specified by the problem
      if(z != 0 && z != zseg+1){
        if(r!= 0 && r!= rseg+1){
          //when r = 1 a phontom point needs to be used as the slope change from
          //0 -> 1 is too large
          //in this case the B field at -1 is assumed to be of the same magnitude as at 1
          //but in the opposite direction
          //this yields the folloeing equation
          if(r==1){
            grids[i + gridStart] += r_aug*(2*grids_old[i+(zseg+2) + gridStart] - 4*grids_old[i + gridStart]) +
              z_aug * (grids_old[i+1+gridStart] - 2*grids_old[i+gridStart] + grids_old[i-1+gridStart]);
          }
          //normal update function
          else{
            grids[i + gridStart] += r_aug*((1+(1/(2*i)))*grids_old[i+(zseg+2) + gridStart] + (-2-(1/(i*i)))*grids_old[i + gridStart] + (1-(1/(2*i)))*grids_old[i-(zseg+2) + gridStart])
              +z_aug*(grids_old[i+1+gridStart] - 2*grids_old[i+gridStart] + grids_old[i-1+gridStart]);
          }
        }
      }
    }
    steps++;
    //wait for all threads
    __syncthreads();
  }
}

int main(){
  //string label in input FILE
  //no use just for clarity
  char label[256];

  double *Il_h, *Ir_h, *rlength_h, *eta_h, *tstep_h, *ldr_h, *ldz_h, *zlength_h, *dI_h, *r_aug_h, *z_aug_h, *grids_h;
  int *rseg_h, *zseg_h;
  long int *totsteps_h, *grid_size_h, *offsets_h;

  double *Il_d, *rlength_d, *ldr_d, *dI_d, *r_aug_d, *z_aug_d, *grids_d, *grids_old_d;
  int *rseg_d, *zseg_d;
  long int *totsteps_d, *grid_size_d, *offsets_d;

  int testcases;
  //printf("How many test cases? ");
  //gets number of test cases
  scanf("%d", &testcases);

  //these sizes are used for malloc ops
  //the size is based on the amount of test cases and the
  //type of data this array holds
  size_t doubleSize = testcases*sizeof(double);
  size_t intSize = testcases*sizeof(int);
  size_t longIntSize = testcases*sizeof(long int);

  //host memory allocation
  Il_h = (double*)malloc(doubleSize);
  Ir_h = (double*)malloc(doubleSize);
  rlength_h = (double*)malloc(doubleSize);
  eta_h = (double*)malloc(doubleSize);
  tstep_h = (double*)malloc(doubleSize);
  ldr_h = (double*)malloc(doubleSize);
  ldz_h = (double*)malloc(doubleSize);
  zlength_h = (double*)malloc(doubleSize);
  dI_h = (double*)malloc(doubleSize);
  r_aug_h = (double*)malloc(doubleSize);
  z_aug_h = (double*)malloc(doubleSize);
  rseg_h = (int*)malloc(intSize);
  zseg_h = (int*)malloc(intSize);
  totsteps_h = (long int*)malloc(longIntSize);
  grid_size_h = (long int*)malloc(longIntSize);
  offsets_h = (long int*)malloc(longIntSize);

  //device memory allocation
  cudaMalloc(&Il_d, doubleSize);
  cudaMalloc(&rlength_d, doubleSize);
  cudaMalloc(&ldr_d, doubleSize);
  cudaMalloc(&dI_d, doubleSize);
  cudaMalloc(&r_aug_d, doubleSize);
  cudaMalloc(&z_aug_d, doubleSize);
  cudaMalloc(&rseg_d, intSize);
  cudaMalloc(&zseg_d, intSize);
  cudaMalloc(&totsteps_d, longIntSize);
  cudaMalloc(&grid_size_d, longIntSize);
  cudaMalloc(&offsets_d, longIntSize);

  long int total_grid_size = 0;
  for(int counter = 0; counter < testcases; counter++){
    //label is unused
    scanf("%s", &label);
    //value of current when z = 0
    //printf("What is your left I? ");
    scanf("%lf", Il_h + counter);
    //value of current when z = zMax
    //printf("What is your right I? ");
    scanf("%lf", Ir_h + counter);
    //total r value
    //printf("What is the radius of your rod? ");
    scanf("%lf", rlength_h + counter);
    //total z value
    //printf("What is the length of your rod? ");
    scanf("%lf", zlength_h + counter);
    //value of eta (proportional to diffusivity)
    //printf("What is eta? ");
    scanf("%lf", eta_h + counter);
    //number of r segments
    //printf("How many segments would you like per radius? ");
    scanf("%d", rseg_h + counter);
    //number of z segments
    //printf("How many segments would you like per length? ");
    scanf("%d", zseg_h + counter);
    //length of each r segment
    ldr_h[counter] = rlength_h[counter]/(rseg_h[counter]+1);
    //length of each z segment
    ldz_h[counter] = zlength_h[counter]/(zseg_h[counter]+1);
    double smallest = ldr_h[counter];
    if(ldz_h[counter] < smallest)
      smallest = ldz_h[counter];
    //determines tstep that ensures stability
    //0.125 is derived from 2^-(num dimensions + 1)
    tstep_h[counter] = 0.125*smallest*smallest*mu0/eta_h[counter];
    //gets total run time
    double tottime;
    //printf("How long would you like to run? ");
    scanf("%lf", &tottime);
    //total steps is an integer that truncates remainder
    totsteps_h[counter] = tottime/tstep_h[counter];
    //dI is change in I per each z segment
    //zseg + 2 to account for boundary conditions
    dI_h[counter] = (Ir_h[counter] - Il_h[counter]) / (zseg_h[counter]+2);
    //r and z aug are used for updating values for increasing time
    //as per numerical update
    r_aug_h[counter] = eta_h[counter]*tstep_h[counter]/(mu0*ldr_h[counter]*ldr_h[counter]);
    z_aug_h[counter] = eta_h[counter]*tstep_h[counter]/(mu0*ldz_h[counter]*ldz_h[counter]);
    //grid size accounts for boundary conditions
    grid_size_h[counter] = (rseg_h[counter] + 2)*(zseg_h[counter] + 2);
    //offsets says where in memory this problems data starts
    //more info above
    offsets_h[counter] = total_grid_size;
    //Total grid size is the size of all needed grid memory
    total_grid_size += grid_size_h[counter];
  }


  size_t gridsSize = total_grid_size*sizeof(double);

  //allocate variable amount of grid memory on both host and device
  grids_h = (double*)malloc(gridsSize);
  cudaMalloc(&grids_d, gridsSize);
  cudaMalloc(&grids_old_d, gridsSize);

  //copy memory down to device
  cudaMemcpy(Il_d, Il_h, doubleSize, cudaMemcpyHostToDevice);
  cudaMemcpy(rlength_d, rlength_h, doubleSize, cudaMemcpyHostToDevice);
  cudaMemcpy(ldr_d, ldr_h, doubleSize, cudaMemcpyHostToDevice);
  cudaMemcpy(dI_d, dI_h, doubleSize, cudaMemcpyHostToDevice);
  cudaMemcpy(r_aug_d, r_aug_h, doubleSize, cudaMemcpyHostToDevice);
  cudaMemcpy(z_aug_d, z_aug_h, doubleSize, cudaMemcpyHostToDevice);
  cudaMemcpy(rseg_d, rseg_h, intSize, cudaMemcpyHostToDevice);
  cudaMemcpy(zseg_d, zseg_h, intSize, cudaMemcpyHostToDevice);
  cudaMemcpy(totsteps_d, totsteps_h, longIntSize, cudaMemcpyHostToDevice);
  cudaMemcpy(grid_size_d, grid_size_h, longIntSize, cudaMemcpyHostToDevice);
  cudaMemcpy(offsets_d, offsets_h, longIntSize, cudaMemcpyHostToDevice);

  //initialize all problems with initial values as per Dr. Sankaran's instructions
  init<<<testcases,threadsPerBlock>>>(grids_d, Il_d, dI_d, ldr_d, rlength_d, grid_size_d, zseg_d, offsets_d);
  //wait for all device implementation to finish
  cudaDeviceSynchronize();

  //copy initial grid data back to host
  cudaMemcpy(grids_h, grids_d, gridsSize, cudaMemcpyDeviceToHost);


  //outputs all problems initial conditions
  FILE *myfile;
  long int i;
  for(int counter = 0; counter < testcases; counter++){
    long int gridStart = offsets_h[counter];
    char init[] = "initXX.txt";
    int tens = counter / 10;
    int ones = counter % 10;
    init[4] = tens + '0';
    init[5] = ones + '0';
    myfile = fopen(init, "w");
    for(i = 0; i< zseg_h[counter]+1; i++)
      fprintf(myfile, "%lf ", i*ldz_h[counter]);
    fprintf(myfile, "%lf\n", i*ldz_h[counter]);

    for(i = 0; i< rseg_h[counter]+1; i++)
      fprintf(myfile, "%lf ", i*ldr_h[counter]);
    fprintf(myfile, "%lf\n", i*ldr_h[counter]);

    for(i = 0; i< (rseg_h[counter] + 2)*(zseg_h[counter]+2); i++){
      if(i%(zseg_h[counter]+2)==zseg_h[counter]+1)
        fprintf(myfile, "%lf\n", grids_h[i + gridStart]);
      else
        fprintf(myfile, "%lf ", grids_h[i + gridStart]);
    }

    fclose(myfile);
  }



  clock_t begin, end;
  double time_spent;
  begin = clock();

  //run
  //size capped at MaxGridPoints, grid points per problem
  run<<<testcases,threadsPerBlock, MaxGridPoints*sizeof(double)>>>(grids_d, r_aug_d, z_aug_d, totsteps_d, grid_size_d, rseg_d, zseg_d, offsets_d, grids_old_d);
  //wait for all device implementation to finish
  cudaDeviceSynchronize();

  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  //copy final grid data back to host
  cudaMemcpy(grids_h, grids_d, gridsSize, cudaMemcpyDeviceToHost);

  //output all problems final grid data
  for(int counter = 0; counter < testcases; counter++){
    long int gridStart = offsets_h[counter];
    char init[] = "resXX.txt";
    int tens = counter / 10;
    int ones = counter % 10;
    init[3] = tens + '0';
    init[4] = ones + '0';
    myfile = fopen(init, "w");
    for(i = 0; i< zseg_h[counter]+1; i++)
      fprintf(myfile, "%lf ", i*ldz_h[counter]);
    fprintf(myfile, "%lf\n", i*ldz_h[counter]);

    for(i = 0; i< rseg_h[counter]+1; i++)
      fprintf(myfile, "%lf ", i*ldr_h[counter]);
    fprintf(myfile, "%lf\n", i*ldr_h[counter]);

    for(i = 0; i< (rseg_h[counter] + 2)*(zseg_h[counter]+2); i++){
      if(i%(zseg_h[counter]+2)==zseg_h[counter]+1)
        fprintf(myfile, "%lf\n", grids_h[i + gridStart]);
      else
        fprintf(myfile, "%lf ", grids_h[i + gridStart]);
    }

    fclose(myfile);
  }

  //free host memory
  free(Il_h);
  free(Ir_h);
  free(rlength_h);
  free(eta_h);
  free(tstep_h);
  free(ldr_h);
  free(ldz_h);
  free(zlength_h);
  free(dI_h);
  free(r_aug_h);
  free(z_aug_h);
  free(rseg_h);
  free(zseg_h);
  free(totsteps_h);
  free(grid_size_h);
  free(offsets_h);

  //free device memory
  cudaFree(Il_d);
  cudaFree(rlength_d);
  cudaFree(ldr_d);
  cudaFree(dI_d);
  cudaFree(r_aug_d);
  cudaFree(z_aug_d);
  cudaFree(rseg_d);
  cudaFree(zseg_d);
  cudaFree(grids_d);
  cudaFree(grids_old_d);
  cudaFree(totsteps_d);
  cudaFree(grid_size_d);
  cudaFree(offsets_d);



  printf("\n------------------------------------\nExecution took: %lf sec\n", time_spent);
  return 0;
}
