#ifndef _MASSCONSERVATION_
#define _MASSCONSERVATION_

#include "globals.h"

//NOTE there are nr+1 points per row of finite volumes
//also mass is a 2d grid compressed to one dimension
#define mOld(i,k) mOld[k*(nr+1)+i]
#define mNew(i,k) mNew[k*(nr+1)+i]
#define mvR(i,k) mvR[k*(nr+1)+i]
#define mvZ(i,k) mvZ[k*(nr+1)+i]
#define massSource(i,k) massSource[k*(nr+1)+i]
//Shared position is defined relative to the current position
//+1 accounts for the -1 in the declaration of i and k
//this -1 is specific for finite volume
#define mvRShared(xOffset,yOffset) mvR_s[(threadIdx.y+yOffset+1)*(R_EVALS_PER_BLOCK+2) + (threadIdx.x+xOffset+1)]
#define mvZShared(xOffset,yOffset) mvZ_s[(threadIdx.y+yOffset+1)*(R_EVALS_PER_BLOCK+2) + (threadIdx.x+xOffset+1)]

__global__ void updateMass(double *mOld, double*mNew, double *mvR, double *mvZ, double *massSource, int nr, int nz, double dr, double dz, double dt){
  extern __shared__ double mv_s[];
  double *mvR_s, *mvZ_s;
  mvR_s = mv_s;
  mvZ_s = mv_s + ((R_EVALS_PER_BLOCK+2)*(Z_EVALS_PER_BLOCK+2));
  //in this case when i or k = 0 it is not a boundary position rather
  //an internal volume next to the boundary
  //to get a boundary volume we subtract 1 to enable the creation of a -1
  //"ghost point"
  int i = blockIdx.x*R_EVALS_PER_BLOCK + threadIdx.x - 1;
  int k = blockIdx.y*Z_EVALS_PER_BLOCK + threadIdx.y - 1;

//Copy data into shared memory for speedup
  if((i == -1 || i == nr+1)){//r ghost point
    if(k != -1 && k != nz + 1){//not corner point
      mvRShared(0,0) = (i == -1 ? mvR(0,k) : mvR(nr,k));
      mvZShared(0,0) = (i == -1 ? mvZ(0,k) : mvZ(nr,k));
    }
  }

  else if(k == -1 || k == nz+1){//z ghost point
    if(i != -1 && i != nr + 1){//not corner point
      mvRShared(0,0) = (k == -1 ? mvR(i,0) : mvR(i,nz));
      mvZShared(0,0) = (k == -1 ? mvZ(i,0) : mvZ(i,nz));
    }
  }

  else if(i > -1 && i< nr + 1 && k > -1 && k < nz + 1){//inside grid
    mvRShared(0,0) = mvR(i,k);
    mvZShared(0,0) = mvZ(i,k);
  }

  __syncthreads();//Wait for all memory copy to be done

  if(i > -1 && i< nr + 1 && k > -1 && k < nz + 1){//inside gird
    if(threadIdx.x != 0 && threadIdx.x != R_EVALS_PER_BLOCK + 1 && threadIdx.y != 0 && threadIdx.y != Z_EVALS_PER_BLOCK + 1){//not halo points
      //Calculate fluxes
      double mvZBot, mvZTop, mvRLeft, mvRRight;


      //TODO add diffusivity error correction
      mvZBot = 0.5*(mvZShared(0,(-1)) + mvZShared(0,0));
      mvZTop = 0.5*(mvZShared(0,0) + mvZShared(0,1));
      mvRLeft = 0.5*(mvRShared((-1),0) + mvRShared(0,0));
      mvRRight = 0.5*(mvRShared(0,0) + mvRShared(1,0));

      mNew(i,k) = mOld(i,k) - dt*((mvRRight-mvRLeft)/dr + (mvZTop-mvZBot)/dz - massSource(i,k));
    }//end not halo points
  }//end inside grid
}

//TODO test function
void getMass(double *mOldP, double *mNewP, double *mvRP, double *mvZP, double *massSourceP,
  double *mOldN, double * mNewN, double *mvRN, double *mvZN, double *massSourceN,
  int nr, int nz, double dr, double dz, double dt,
  dim3 centerGridWHalosBlockDim, dim3 centerGridWHalosThreadDim, size_t centerGridSize){

      updateMass<<<centerGridWHalosBlockDim,centerGridWHalosThreadDim,2*centerGridSize>>>
      (mOldP, mNewP, mvRP, mvZP, massSourceP, nr, nz, dr, dz, dt);//Update mass positives

      updateMass<<<centerGridWHalosBlockDim,centerGridWHalosThreadDim,2*centerGridSize>>>
      (mOldN, mNewN, mvRN, mvZN, massSourceN, nr, nz, dr, dz, dt);//Update mass negatives
  }


#endif
