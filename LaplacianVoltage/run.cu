#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#define PI 3.14159265359

#define grid(i,k,nr) k*nr+i
#define omega 1.5
#define RelativeError 1e-3
#define epsilon 1e-12
#define nMax 256

#define zEvalsPerBlock 16
#define rEvalsPerBlock 16

double Besseli0(double x){
  //returns modified bessel function I_0 for any real x
  double ax, ans;
  double y;

  if((ax=fabs(x)) < 3.75){
    y = x/3.75;
    y*=y;
    ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
  }
  else {
    y=3.75/ax;
    ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1+y*0.392377e-2))))))));
  }
  return ans;
}

double B(int n){
  return 200.0/(n*PI*Besseli0(n*PI/2));
}

double V(double r, double z, double lz){
  double sum = 0;
  for (int n = 1; n <= nMax; n++) {
    double b = B(n);
    sum += b*Besseli0(n*PI/lz*r)*sin(n*PI/lz*z);
  }
  return sum;
}

//-----------------------------------------------------------------------------
__global__ void update(double *volt_old, double * volt_new, bool isRed, bool *converge, int nr, int nz, double dr, double dz, int blockz){
  extern __shared__ double volt_s[];
  int i = blockIdx.x * zEvalsPerBlock + threadIdx.x;//x position index
  int k = blockIdx.y * rEvalsPerBlock + threadIdx.y;//y position index
  int sharedPos = threadIdx.y * (rEvalsPerBlock+2) + threadIdx.x;//Making 2D into 1D for the shared memory position in this block
  int blockPos = 2*(blockIdx.x * blockz + blockIdx.y);
  if(isRed){//Could have just as well been !isRed
    blockPos++;//Because each block has two convergence flags, need to only update one of the two
  }
  if(i< nr+2 && k < nz+2 && i != 0){//Within the domain of the grid
    volt_s[sharedPos] = volt_old[grid(i,k,(nr+2))];
  }
  //Because value of center of nozzle is a floating potential copy the value that is above it
  if(i==0){
    volt_s[sharedPos] = volt_old[grid(1,k,(nr+2))];//Zero gradient between r=0 and r=dr
  }
  __syncthreads();//This is to ensure that all the threads have copied values from the previous iteration to shared memory
  if((i%2 == k%2) == isRed){//Red or not Red?
    converge[blockPos] = true;//Default. Then all you need is one 'false' to force another iteration
  }
  __syncthreads();//ensures converge has been set to true for all threads

  if((i%2 == k%2) == isRed){//Red or not Red?
    if(threadIdx.x != 0 && threadIdx.y != 0 && threadIdx.x != zEvalsPerBlock + 1 && threadIdx.y != rEvalsPerBlock + 1){//not halo points
      if(i != 0 && i < nr+1 && k != 0 && k < nz+1){//not boundaries

        volt_new[grid(i,k,(nr+2))] = (1-omega)*volt_s[sharedPos];//copy a weighted fraction of the old
        //Then update with the remaining fraction with the new
        volt_new[grid(i,k,(nr+2))] += omega *(
          volt_s[sharedPos-(rEvalsPerBlock+2)]*dr*dr/(2*(dr*dr+dz*dz)) + //Bottom
          volt_s[sharedPos+rEvalsPerBlock+2]*dr*dr/(2*(dr*dr+dz*dz)) + //Top
          volt_s[sharedPos-1]*dz*dz/(2*(dr*dr+dz*dz))*(1-(1/(2*i))) + //Left
          volt_s[sharedPos+1]*dz*dz/(2*(dr*dr+dz*dz))*(1+(1/(2*i))) //Right
        );

        //Convergence check
        double relChange = fabs((volt_new[grid(i,k,(nr+2))] - volt_s[sharedPos])/(volt_s[sharedPos] + epsilon));
        if(relChange > RelativeError/max(nr,nz)){
          converge[blockPos] = false;
        }//end converge check

      }//end of not boundaries
    }//end of not halo
  }//end of red/black
}//end of update


__global__ void getEFields(double *Er, double *Ez, double *volt, double dr, double dz, int nr){
  extern __shared__ double volt_s[];
  int i = blockIdx.x * zEvalsPerBlock + threadIdx.x;//x position index
  int k = blockIdx.y * rEvalsPerBlock + threadIdx.y;//y position index
  int sharedPos = threadIdx.y * (rEvalsPerBlock+1) + threadIdx.x;//Making 2D into 1D for the shared memory position in this block
  volt_s[sharedPos] = volt[grid(i,k,(nr+2))];//Bottom Left
  volt_s[sharedPos + 1] = volt[grid((i+1),k,(nr+2))];//Bottom right
  volt_s[sharedPos + rEvalsPerBlock + 1] = volt[grid(i,k,(nr+2))];//Top Left
  volt_s[sharedPos + rEvalsPerBlock + 1 + 1] = volt[grid(i,k,(nr+2))];//Top right
  double vtop, vbot, vleft, vright;
  vtop = (volt_s[sharedPos + rEvalsPerBlock + 1] + volt_s[sharedPos + rEvalsPerBlock + 1 + 1])/2;
  vbot = (volt_s[sharedPos] + volt_s[sharedPos + 1])/2;
  vleft = (volt_s[sharedPos] + volt_s[sharedPos + rEvalsPerBlock + 1])/2;
  vright = (volt_s[sharedPos + 1] + volt_s[sharedPos + rEvalsPerBlock + 1 + 1])/2;
  Er[grid(i,k,(nr+1))] = -(vright-vleft)/dr;
  Ez[grid(i,k,(nr+1))] = -(vtop-vbot)/dz;
}

//-----------------------------------------------------------------------------
int main(){

  /*
    NOTE: i index refers to radial direction and it is treated equivalent to x
          k index refers to axial direction
          left and right in code refer to change in r while top and bottom refers
          to change in z
          This changes in plotting where z is plotted on the horizontal axis and
          r is plotted on the vertical axis
  */

  //Get boundary conditions
  double vleft, vright;
  scanf("%*s");
  scanf("%lf", &vleft);
  scanf("%*s");
  scanf("%lf", &vright);

  //Get segmentation info
  int nz, nr;
  scanf("%*s");
  scanf("%d", &nz);
  scanf("%*s");
  scanf("%d", &nr);

  //Get length info
  double lz, lr;
  scanf("%*s");
  scanf("%lf", &lz);
  scanf("%*s");
  scanf("%lf", &lr);

  //Calculate segment length
  double dz = lz/(nz+1);
  double dr = lr/(nr+1);

  //Allocate grid memory on both host and device
  double *volt_h; //for host
  double *volt_d_old, *volt_d_new; //for device
  double *error;
  size_t gridSize = (nr+2)*(nz+2)*sizeof(double);
  cudaMalloc(&volt_d_new, gridSize);
  cudaMalloc(&volt_d_old, gridSize);
  volt_h = (double*)malloc(gridSize);
  error = (double*)malloc(gridSize);

//Boundary conditions
  //get change in voltage per change in z
  double VdeltaZ = 100.0/(nz+1);

//Left and right
  for(int i = 0; i< nr + 2; i++){
    volt_h[grid(i,0,(nr+2))] = vleft;
    volt_h[grid(i,(nz+1),(nr+2))] = vright;
  }

//Top
  for(int k=0; k< nz+2; k++){
    volt_h[grid((nr+1),k,(nr+2))] = 100 - VdeltaZ*k;
  }

//initial guess
for(int i = 1; i< nr + 1; i++){
  for(int k = 1; k< nz + 1; k++){
    volt_h[grid(i,k,(nr+2))] = (vleft + vright)/2;
  }
}


  //copy memory down to devicce
  cudaMemcpy(volt_d_new, volt_h, gridSize, cudaMemcpyHostToDevice);


  //Get block dimenions
  dim3 threadDim(zEvalsPerBlock + 2, rEvalsPerBlock + 2);// +2 accounts for halo points
  int blockr = 1 + ((nr-1)/rEvalsPerBlock);//nr is the number of interior r points
  int blockz = 1 + ((nz-1)/zEvalsPerBlock);//nz is the number of interior z points
  dim3 blockDim(blockr, blockz);

  //Allocate converge check memory
  bool *converge_d;
  bool *converge_h;
  size_t convSize = 2*blockr*blockz*sizeof(bool);//Each block has a red conv pos and a black conv pos therefore 2*numBlocks
  converge_h = (bool*)malloc(convSize);
  cudaMalloc(&converge_d, convSize);

  int steps = 0;
  bool didConverge = false;
  printf("Blocks: %d\n", blockr*blockz);
  while(!didConverge){
    cudaMemcpy(volt_d_old, volt_d_new, gridSize, cudaMemcpyDeviceToDevice);
    //Evaluate red blocks
    update<<<blockDim,threadDim,(zEvalsPerBlock+2)*(rEvalsPerBlock+2)*sizeof(double)>>>(volt_d_old, volt_d_new, true, converge_d, nr, nz, dr, dz, blockz);
    cudaMemcpy(volt_d_old, volt_d_new, gridSize, cudaMemcpyDeviceToDevice);
    //Evaluate black blocks
    update<<<blockDim,threadDim,(zEvalsPerBlock+2)*(rEvalsPerBlock+2)*sizeof(double)>>>(volt_d_old, volt_d_new, false, converge_d, nr, nz, dr, dz, blockz);
    //copy back converge check
    cudaMemcpy(converge_h, converge_d, convSize, cudaMemcpyDeviceToHost);
    steps++;

    //all converge must be true
    didConverge = converge_h[0];
    for(int i = 1; i< 2*blockr*blockz; i++){
      didConverge = didConverge && converge_h[i];
    }
  }//converged
  printf("\nconverged in %d steps.\n", steps);
  cudaMemcpy(volt_h, volt_d_new, gridSize, cudaMemcpyDeviceToHost);

  for(int k = 0; k< nz+2; k++){
    volt_h[grid(0,k,(nr+2))] = volt_h[grid(1,k,(nr+2))];
  }
//---------------------------------------------------------------------------//
//This section deals with E fields

/*
  size_t Esize = (nr+1)*(nz+1)*sizeof(double);
  double *Er_d, *Ez_d, *Er_h, *Ez_h;
  cudaMalloc(&Er_d, Esize);
  cudaMalloc(&Ez_d, Esize);
  Er_h = (double*)malloc(Esize);
  Ez_h = (double*)malloc(Esize);

  //Get block dimenions
  dim3 enumthreadDim(zEvalsPerBlock, rEvalsPerBlock);
  int eblockr = 1 + ((nr+1-1)/rEvalsPerBlock);//nr is the number of interior r points
  int eblockz = 1 + ((nz+1-1)/zEvalsPerBlock);//nz is the number of interior z points
  dim3 eblockDim(eblockr, eblockz);

  getEFields<<<eblockDim,enumthreadDim,(rEvalsPerBlock + 1)*(zEvalsPerBlock + 1)*sizeof(double)>>>();
  cudaMemcpy(Er_h, Er_d, Esize, cudaMemcpyDeviceToHost);
  cudaMemcpy(Ez_h, Ez_d, Esize, cudaMemcpyDeviceToHost);
*/


//----------------------------------------------------------------------------//
  //Output data for contour maps plotly
  FILE *myfile, *err;
  myfile = fopen("voltageRes.txt", "w");
  err = fopen("voltageErr.txt", "w");
  int i, j;
  for(i = 0; i< nz+1; i++){
    fprintf(myfile, "%lf ", i*dz);
    fprintf(err, "%lf ", i*dz);
  }
  fprintf(myfile, "%lf\n", i*dz);//This is the nz+1 point.
  fprintf(err, "%lf\n", i*dz);//This is the nz+1 point.
  for(j = 0; j<nr+1; j++){
    fprintf(myfile, "%lf ", j*dr);
    fprintf(err, "%lf ", j*dr);
  }
  fprintf(myfile, "%lf\n", j*dr);//This is the nr+1 point.
  fprintf(err, "%lf\n", j*dr);//This is the nr+1 point.


  for(j = 0; j<nr+2; j++){
    for(i = 0; i< nz+1; i++){
      if(i==0 || j==nr+1){
        error[grid(j,i,(nr+2))] = 0;
      }
      else{
        error[grid(j,i,(nr+2))] =volt_h[grid(j,i,(nr+2))] - V(j*dr, i*dz, lz);
      }
      fprintf(myfile, "%lf ", volt_h[grid(j,i,(nr+2))]);
      fprintf(err, "%lf ", error[grid(j,i,(nr+2))]);
    }
    error[grid(j,i,(nr+2))] = 0;
    fprintf(myfile, "%lf\n", volt_h[grid(j,i,(nr+2))]);//nz+1 point for each j
    fprintf(err, "%lf\n", error[grid(j,i,(nr+2))]);//nz+1 point for each j
  }
  fprintf(myfile, "%lf %lf\n", dz, dr);
  fclose(myfile);
  fclose(err);

  cudaFree(volt_d_new);
  cudaFree(volt_d_old);
  cudaFree(converge_d);
  //cudaFree(Er_d);
  //cudaFree(Ez_d);
  //free(Er_h);
  //free(Ez_h);
  free(volt_h);
  free(converge_h);
}
