#ifndef _MOMENTUMCONSERVATION_H_
#define _MOMENTUMCONSERVATION_H_
//TODO ghost cells
#include "globals.h"

//NOTE there are nr+1 points per row of finite volumes
//also these are 2d grids compressed to one dimension
#define massOld(i,k) mOld[k*(nr+1)+i]
#define mvROld(i,k) mvROld[k*(nr+1)+i]
#define mvZOld(i,k) mvZOld[k*(nr+1)+i]
#define mvRNew(i,k) mvRNew[k*(nr+1)+i]
#define mvZNew(i,k) mvZNew[k*(nr+1)+i]
#define mvSource(i,k) mvSource[k*(nr+1)+i]
//Shared position is defined relative to the current position
//+1 accounts for the -1 in the declaration of i and k
//this -1 is specific for finite volume
#define mvROldShared(xOffset,yOffset) mvROld_s[(threadIdx.y+yOffset+1)*(R_EVALS_PER_BLOCK+2) + (threadIdx.x+xOffset+1)]
#define mvZOldShared(xOffset,yOffset) mvZOld_s[(threadIdx.y+yOffset+1)*(R_EVALS_PER_BLOCK+2) + (threadIdx.x+xOffset+1)]
#define mOldShared(xOffset,yOffset) mass_s[(threadIdx.y+yOffset+1)*(R_EVALS_PER_BLOCK+2) + (threadIdx.x+xOffset+1)]

__global__ void updateMvRhat(double *mOld, double*mvRNew, double *mvROld, double *mvZOld, double *mvSource, int nr, int nz, double dr, double dz, double dt, bool polarity, int atomicMass){
  extern __shared__ double U_s[];
  double *mvROld_s, *mvZOld_s, *mass_s;
  mvROld_s = U_s;
  mvZOld_s = U_s + ((R_EVALS_PER_BLOCK+2)*(Z_EVALS_PER_BLOCK+2));
  mass_s = U_s + (2*(R_EVALS_PER_BLOCK+2)*(Z_EVALS_PER_BLOCK+2));
  //in this case when i or k = 0 it is not a boundary position rather
  //an internal volume next to the boundary
  //to get a boundary volume we subtract 1 to enable the creation of a -1
  //"ghost point"
  int i = blockIdx.x*R_EVALS_PER_BLOCK + threadIdx.x - 1;
  int k = blockIdx.y*Z_EVALS_PER_BLOCK + threadIdx.y - 1;

  double Mparticle;
  double Temperature;
  if(polarity){
    Mparticle = atomicMass*AMU;
    Temperature = Tp;
  }
  else{
    Mparticle = atomicMass*AMU;
    Temperature = Tn;
  }

//Copy data into shared memory for speedup
  if((i == -1 || i == nr+1)){//r ghost point
    if(k != -1 && k != nz + 1){//not corner point
      mvROldShared(0,0) = (i == -1 ? mvROld(0,k) : mvROld(nr,k));
      mvZOldShared(0,0) = (i == -1 ? mvZOld(0,k) : mvZOld(nr,k));
      mOldShared(0,0) = (i == -1 ? massOld(0,k) : massOld(nr,k));
    }
  }

  else if(k == -1 || k == nz+1){//z ghost point
    if(i != -1 && i != nr + 1){//not corner point
      mvROldShared(0,0) = (k == -1 ? mvROld(i,0) : mvROld(i,nz));
      mvZOldShared(0,0) = (k == -1 ? mvZOld(i,0) : mvZOld(i,nz));
      mOldShared(0,0) = (k == -1 ? massOld(i,0) : massOld(i,nz));
    }
  }

  else if(i > -1 && i< nr + 1 && k > -1 && k < nz + 1){//inside grid
    mvROldShared(0,0) = mvROld(i,k);
    mvZOldShared(0,0) = mvZOld(i,k);
    mOldShared(0,0) = massOld(i,k);
  }

  __syncthreads();//Wait for all memory copy to be done

  if(i > -1 && i< nr + 1 && k > -1 && k < nz + 1){//inside gird
    if(threadIdx.x != 0 && threadIdx.x != R_EVALS_PER_BLOCK + 1 && threadIdx.y != 0 && threadIdx.y != Z_EVALS_PER_BLOCK + 1){//not halo points
      //Calculate fluxes
      double botFace, topFace, leftFace, rightFace;


      //TODO add diffusivity error correction
      leftFace = 0.5*((mvROldShared((-1),0)*mvROldShared((-1),0)/mOldShared((-1),0)+mOldShared((-1),0)*kBoltz*Temperature/Mparticle)+(mvROldShared(0,0)*mvROldShared(0,0)/mOldShared(0,0)+mOldShared(0,0)*kBoltz*Temperature/Mparticle));
      rightFace = 0.5*((mvROldShared(0,0)*mvROldShared(0,0)/mOldShared(0,0)+mOldShared(0,0)*kBoltz*Temperature/Mparticle)+(mvROldShared(1,0)*mvROldShared(1,0)/mOldShared(1,0)+mOldShared(1,0)*kBoltz*Temperature/Mparticle));
      botFace = 0.5*((mvROldShared(0,(-1))*mvZOldShared(0,(-1))/mOldShared(0,(-1)))+(mvROldShared(0,0)*mvZOldShared(0,0)/mOldShared(0,0)));
      topFace = 0.5*((mvROldShared(0,0)*mvZOldShared(0,0)/mOldShared(0,0))+(mvROldShared(0,1)*mvZOldShared(0,1)/mOldShared(0,1)));



      mvRNew(i,k) = mvROld(i,k) - dt*((rightFace-leftFace)/dr + (topFace-botFace)/dz - mvSource(i,k));
    }//end not halo points
  }//end inside grid
}

__global__ void updateMvZhat(double *mOld, double*mvZNew, double *mvROld, double *mvZOld, double *mvSource, int nr, int nz, double dr, double dz, double dt, bool polarity, int atomicMass){
  extern __shared__ double U_s[];
  double *mvROld_s, *mvZOld_s, *mass_s;
  mvROld_s = U_s;
  mvZOld_s = U_s + ((R_EVALS_PER_BLOCK+2)*(Z_EVALS_PER_BLOCK+2));
  mass_s = U_s + (2*(R_EVALS_PER_BLOCK+2)*(Z_EVALS_PER_BLOCK+2));
  //in this case when i or k = 0 it is not a boundary position rather
  //an internal volume next to the boundary
  //to get a boundary volume we subtract 1 to enable the creation of a -1
  //"ghost point"
  int i = blockIdx.x*R_EVALS_PER_BLOCK + threadIdx.x - 1;
  int k = blockIdx.y*Z_EVALS_PER_BLOCK + threadIdx.y - 1;

  double Mparticle;
  double Temperature;
  if(polarity){
    Mparticle = atomicMass*AMU;
    Temperature = Tp;
  }
  else{
    Mparticle = atomicMass*AMU;
    Temperature = Tn;
  }

//Copy data into shared memory for speedup
  if((i == -1 || i == nr+1)){//r ghost point
    if(k != -1 && k != nz + 1){//not corner point
      mvROldShared(0,0) = (i == -1 ? mvROld(0,k) : mvROld(nr,k));
      mvZOldShared(0,0) = (i == -1 ? mvZOld(0,k) : mvZOld(nr,k));
      mOldShared(0,0) = (i == -1 ? massOld(0,k) : massOld(nr,k));
    }
  }

  else if(k == -1 || k == nz+1){//z ghost point
    if(i != -1 && i != nr + 1){//not corner point
      mvROldShared(0,0) = (k == -1 ? mvROld(i,0) : mvROld(i,nz));
      mvZOldShared(0,0) = (k == -1 ? mvZOld(i,0) : mvZOld(i,nz));
      mOldShared(0,0) = (k == -1 ? massOld(i,0) : massOld(i,nz));
    }
  }

  else if(i > -1 && i< nr + 1 && k > -1 && k < nz + 1){//inside grid
    mvROldShared(0,0) = mvROld(i,k);
    mvZOldShared(0,0) = mvZOld(i,k);
    mOldShared(0,0) = massOld(i,k);
  }

  __syncthreads();//Wait for all memory copy to be done

  if(i > -1 && i< nr + 1 && k > -1 && k < nz + 1){//inside gird
    if(threadIdx.x != 0 && threadIdx.x != R_EVALS_PER_BLOCK + 1 && threadIdx.y != 0 && threadIdx.y != Z_EVALS_PER_BLOCK + 1){//not halo points
      //Calculate fluxes
      double botFace, topFace, leftFace, rightFace;


      //TODO add diffusivity error correction
      leftFace = 0.5*((mvZOldShared((-1),0)*mvROldShared((-1),0)/mOldShared((-1),0))+(mvZOldShared(0,0)*mvROldShared(0,0)/mOldShared(0,0)));
      rightFace = 0.5*((mvZOldShared(0,0)*mvROldShared(0,0)/mOldShared(0,0))+(mvZOldShared(1,0)*mvROldShared(1,0)/mOldShared(1,0)));
      botFace = 0.5*((mvZOldShared(0,(-1))*mvZOldShared(0,(-1))/mOldShared(0,(-1))+mOldShared(0,(-1))*kBoltz*Temperature/Mparticle)+(mvZOldShared(0,0)*mvZOldShared(0,0)/mOldShared(0,0)+mOldShared(0,0)*kBoltz*Temperature/Mparticle));
      topFace = 0.5*((mvZOldShared(0,0)*mvZOldShared(0,0)/mOldShared(0,0)+mOldShared(0,0)*kBoltz*Temperature/Mparticle)+(mvZOldShared(0,1)*mvZOldShared(0,1)/mOldShared(0,1)+mOldShared(0,1)*kBoltz*Temperature/Mparticle));



      mvZNew(i,k) = mvZOld(i,k) - dt*((rightFace-leftFace)/dr + (topFace-botFace)/dz - mvSource(i,k));
    }//end not halo points
  }//end inside grid
}

//TODO test function
void getMomentum(double *mOldP, double *mvRNewP, double *mvZNewP, double *mvROldP, double *mvZOldP, double *mvRSourceP, double *mvZSourceP,
  double *mOldN, double *mvRNewN, double *mvZNewN, double *mvROldN, double *mvZOldN, double *mvRSourceN, double *mvZSourceN,
  int nr, int nz, double dr, double dz, double dt, int atomicMass,
  dim3 centerGridWHalosBlockDim, dim3 centerGridWHalosThreadDim){

      updateMvRhat<<<centerGridWHalosBlockDim,centerGridWHalosThreadDim,3*(R_EVALS_PER_BLOCK+2)*(Z_EVALS_PER_BLOCK+2)*sizeof(double)>>>
      (mOldP, mvRNewP, mvROldP, mvZOldP, mvRSourceP, nr, nz, dr, dz, dt, true, atomicMass);//Update r hat momentum positives

      updateMvRhat<<<centerGridWHalosBlockDim,centerGridWHalosThreadDim,3*(R_EVALS_PER_BLOCK+2)*(Z_EVALS_PER_BLOCK+2)*sizeof(double)>>>
      (mOldN, mvRNewN, mvROldN, mvZOldN, mvRSourceN, nr, nz, dr, dz, dt, false, atomicMass);//Update r hat momentum negatives

      updateMvZhat<<<centerGridWHalosBlockDim,centerGridWHalosThreadDim,3*(R_EVALS_PER_BLOCK+2)*(Z_EVALS_PER_BLOCK+2)*sizeof(double)>>>
      (mOldP, mvZNewP, mvROldP, mvZOldP, mvZSourceP, nr, nz, dr, dz, dt, true, atomicMass);//Update z hat momentum positives

      updateMvZhat<<<centerGridWHalosBlockDim,centerGridWHalosThreadDim,3*(R_EVALS_PER_BLOCK+2)*(Z_EVALS_PER_BLOCK+2)*sizeof(double)>>>
      (mOldN, mvZNewN, mvROldN, mvZOldN, mvZSourceN, nr, nz, dr, dz, dt, false, atomicMass);//Update z hat momentum negatives

  }


#endif
