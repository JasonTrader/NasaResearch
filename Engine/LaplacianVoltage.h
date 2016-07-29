#ifndef _LAPLACIANVOLTAGE_H_
#define _LAPLACIANVOLTAGE_H_

#include "globals.h"

#define vOld(i,k) voltOld[k*(nr+2)+i]
#define vNew(i,k) voltNew[k*(nr+2)+i]
#define vShared(xOffset,yOffset) volt_s[(threadIdx.y+yOffset)*(R_EVALS_PER_BLOCK+2) + (threadIdx.x+xOffset)]
#define OMEGA 1.5
#define RELATIVE_ERROR 1e-3


//-----------------------------------------------------------------------------
__global__ void updateVoltage(double *voltOld, double * voltNew, bool isRed, bool *converge, int nr, int nz, double dr, double dz, int blockZ){
  extern __shared__ double volt_s[];
  int i = blockIdx.x * Z_EVALS_PER_BLOCK + threadIdx.x;//x position index
  int k = blockIdx.y * R_EVALS_PER_BLOCK + threadIdx.y;//y position index
  int blockPos = 2*(blockIdx.x * blockZ + blockIdx.y);
  if(isRed){//Could have just as well been !isRed
    blockPos++;//Because each block has two convergence flags, need to only update one of the two
  }
  if(i< nr+2 && k < nz+2 && i != 0){//Within the domain of the grid
    vShared(0,0) = vOld(i,k);
  }
  //Because value of center of nozzle is a floating potential copy the value that is above it
  if(i==0){
    vShared(0,0) = vOld(i,k);//Zero gradient between r=0 and r=dr
  }
  __syncthreads();//This is to ensure that all the threads have copied values from the previous iteration to shared memory
  if((i%2 == k%2) == isRed){//Red or not Red?
    converge[blockPos] = true;//Default. Then all you need is one 'false' to force another iteration
  }
  __syncthreads();//ensures converge has been set to true for all threads

  if((i%2 == k%2) == isRed){//Red or not Red?
    if(threadIdx.x != 0 && threadIdx.y != 0 && threadIdx.x != Z_EVALS_PER_BLOCK + 1 && threadIdx.y != R_EVALS_PER_BLOCK + 1){//not halo points
      if(i != 0 && i < nr+1 && k != 0 && k < nz+1){//not boundaries

        vNew(i,k) = (1-OMEGA)*vShared(0,0);//copy a weighted fraction of the old
        //Then update with the remaining fraction with the new
        vNew(i,k) += OMEGA *(
          vShared(0,-1)*dr*dr/(2*(dr*dr+dz*dz)) + //Bottom
          vShared(0,1)*dr*dr/(2*(dr*dr+dz*dz)) + //Top
          vShared(-1,0)*dz*dz/(2*(dr*dr+dz*dz))*(1-(1/(2*i))) + //Left
          vShared(1,0)*dz*dz/(2*(dr*dr+dz*dz))*(1+(1/(2*i))) //Right
        );//TODO implement source sink term

        //Convergence check
        double relChange = fabs((vNew(i,k) - vShared(0,0))/(vShared(0,0) + EPSILON));
        if(relChange > RELATIVE_ERROR/max(nr,nz)){
          converge[blockPos] = false;
        }//end converge check

      }//end of not boundaries
    }//end of not halo
  }//end of red/black
}//end of update

void getNewVoltage(size_t cornerGridSize, size_t convSize, double *voltOld_d, double *voltNew_d,
  double *volt_h, dim3 cornerGridWHalosBlockDim, dim3 cornerGridWHalosThreadDim,
   bool *converge_d, bool *converge_h, int nr, int nz, double dr, double dz,
   int blockR, int blockZ){


  bool didConverge = false;
  while(!didConverge){
    cudaMemcpy(voltOld_d, voltNew_d, cornerGridSize, cudaMemcpyDeviceToDevice);
    //Evaluate red blocks
    updateVoltage<<<cornerGridWHalosBlockDim,cornerGridWHalosThreadDim,(Z_EVALS_PER_BLOCK+2)*(R_EVALS_PER_BLOCK+2)*sizeof(double)>>>(voltOld_d, voltNew_d, true, converge_d, nr, nz, dr, dz, blockZ);
    cudaMemcpy(voltOld_d, voltNew_d, cornerGridSize, cudaMemcpyDeviceToDevice);
    //Evaluate black blocks
    updateVoltage<<<cornerGridWHalosBlockDim,cornerGridWHalosThreadDim,(Z_EVALS_PER_BLOCK+2)*(R_EVALS_PER_BLOCK+2)*sizeof(double)>>>(voltOld_d, voltNew_d, false, converge_d, nr, nz, dr, dz, blockZ);
    //copy back converge check
    cudaMemcpy(converge_h, converge_d, convSize, cudaMemcpyDeviceToHost);

    //all converge must be true
    didConverge = converge_h[0];
    for(int i = 1; i< 2*blockR*blockZ; i++){
      didConverge = didConverge && converge_h[i];
    }
  }//converged

  cudaMemcpy(volt_h, voltNew_d, cornerGridSize, cudaMemcpyDeviceToHost);//Copy results back
}

#endif
