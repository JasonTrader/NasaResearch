#ifndef _MOMENTUMCONSERVATION_H_
#define _MOMENTUMCONSERVATION_H_

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

__global__ void updateMvRhat(double *mOld, double*mvRNew, double *mvROld, double *mvZOld,
  double *mvSource, int nr, int nz, double dr, double dz, double dt, bool polarity,
  int atomicMass, double rin, double rout, double propellantFlowRate, double a){

  extern __shared__ double U_s[];
  //divide up shared memory
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
  if(polarity){//positive ions
    Mparticle = MASS_POSITIVE_ION;
    Temperature = Tp;
  }
  else{//negative ions
    Mparticle = MASS_NEGATIVE_ION;
    Temperature = Tn;
  }

//Copy data into shared memory for speedup
  if((i != -1 && i != nr+1)){//not r ghost point
    if(k == -1 || k == nz+1){//z ghost point
      if(k==-1){//thruster inlet
        mvROldShared(0,0)=0;//no radial velocity at inlet
        mvZOldShared(0,0)=INLET_MOMENTUM((i+0.5));//Based on propellant flow rate
        mOldShared(0,0)=mvZOldShared(0,0)/a;//Incoming mass
      }
      else{//thruster exit
        //use continuous gradient to approximate
        mvROldShared(0,0)=2*mvROld(i,(k-1))-mvROld(i,(k-2));
        mvZOldShared(0,0) = 2*mvZOld(i,(k-1))-mvZOld(i,(k-2));
        mOldShared(0,0) = 2*massOld(i,(k-1))-massOld(i,(k-2));
      }
    }
    else if(i > -1 && i< nr + 1 && k > -1 && k < nz + 1){//inside grid
      //copy to shared memory
      mvROldShared(0,0) = mvROld(i,k);
      mvZOldShared(0,0) = mvZOld(i,k);
      mOldShared(0,0) = massOld(i,k);
    }
  }

  __syncthreads();//Wait for all memory copy to be done

  if(i > -1 && i< nr + 1 && k > -1 && k < nz + 1){//inside gird
    if(threadIdx.x != 0 && threadIdx.x != R_EVALS_PER_BLOCK + 1 && threadIdx.y != 0 && threadIdx.y != Z_EVALS_PER_BLOCK + 1){//not halo points
      //Calculate fluxes
      double botFace, topFace, leftFace, rightFace;
      double Dbot, Dtop, Dleft, Dright;


      if(i==0){//thruster inlet
        if(rin==0){//center
          leftFace = 0;//no streams across center
          rightFace = 0.5*((mvROldShared(0,0)*mvROldShared(0,0)/mOldShared(0,0)+mOldShared(0,0)*kBoltz*Temperature/Mparticle)+(mvROldShared(1,0)*mvROldShared(1,0)/mOldShared(1,0)+mOldShared(1,0)*kBoltz*Temperature/Mparticle));
        }
        else{//center of pipe
          //r is not 0 because center of pipe is insulated
          leftFace=mOldShared(0,0)*kBoltz*Temperature/Mparticle;//'r' is already multiplied in the mOldShared
          rightFace = 0.5*((mvROldShared(0,0)*mvROldShared(0,0)/mOldShared(0,0)+mOldShared(0,0)*kBoltz*Temperature/Mparticle)+(mvROldShared(1,0)*mvROldShared(1,0)/mOldShared(1,0)+mOldShared(1,0)*kBoltz*Temperature/Mparticle));
        }
        //Numerical dampening
        Dleft = 0;
        Dright = 0.5*0.5*((mvROldShared(1,0)/mOldShared((i+1),k)+a)+(mvROldShared(0,0)/mOldShared(i,k)+a))*(mvROldShared((i+1),k)-mvROldShared(i,k));
      }
      else if(i==nr){//thruster top
        leftFace = 0.5*((mvROldShared((-1),0)*mvROldShared((-1),0)/mOldShared((-1),0)+mOldShared((-1),0)*kBoltz*Temperature/Mparticle)+(mvROldShared(0,0)*mvROldShared(0,0)/mOldShared(0,0)+mOldShared(0,0)*kBoltz*Temperature/Mparticle));
        //no streams through top of rocket
        rightFace=mOldShared(0,0)*kBoltz*Temperature/Mparticle;//'r' is already multiplied in the mOldShared
        //Numerical dampening
        Dleft = 0.5*0.5*((mvROldShared(0,0)/mOldShared(i,k)+a)+(mvROldShared((-1),0)/mOldShared((i-1),k)+a))*(mvROldShared(i,k)-mvROldShared((i-1),k));
        Dright = 0;
      }
      else{//inside thruster
        leftFace = 0.5*((mvROldShared((-1),0)*mvROldShared((-1),0)/mOldShared((-1),0)+mOldShared((-1),0)*kBoltz*Temperature/Mparticle)+(mvROldShared(0,0)*mvROldShared(0,0)/mOldShared(0,0)+mOldShared(0,0)*kBoltz*Temperature/Mparticle));
        rightFace = 0.5*((mvROldShared(0,0)*mvROldShared(0,0)/mOldShared(0,0)+mOldShared(0,0)*kBoltz*Temperature/Mparticle)+(mvROldShared(1,0)*mvROldShared(1,0)/mOldShared(1,0)+mOldShared(1,0)*kBoltz*Temperature/Mparticle));
        //Numerical dampening
        Dleft = 0.5*0.5*((mvROldShared(0,0)/mOldShared(i,k)+a)+(mvROldShared((-1),0)/mOldShared((i-1),k)+a))*(mvROldShared(i,k)-mvROldShared((i-1),k));
        Dright = 0.5*0.5*((mvROldShared(1,0)/mOldShared((i+1),k)+a)+(mvROldShared(0,0)/mOldShared(i,k)+a))*(mvROldShared((i+1),k)-mvROldShared(i,k));
      }
      botFace = 0.5*((mvROldShared(0,(-1))*mvZOldShared(0,(-1))/mOldShared(0,(-1)))+(mvROldShared(0,0)*mvZOldShared(0,0)/mOldShared(0,0)));
      topFace = 0.5*((mvROldShared(0,0)*mvZOldShared(0,0)/mOldShared(0,0))+(mvROldShared(0,1)*mvZOldShared(0,1)/mOldShared(0,1)));
      //Numerical dampening
      Dbot = 0.5*0.5*((mvZOldShared(0,0)/mOldShared(i,k)+a)+(mvZOldShared(0,(-1))/mOldShared(i,(k-1))+a))*(mvROldShared(i,k)-mvROldShared(i,(k-1)));
      Dtop = 0.5*0.5*((mvZOldShared(0,1)/mOldShared(i,(k+1))+a)+(mvZOldShared(0,0)/mOldShared(i,k)+a))*(mvROldShared(i,(k+1))-mvROldShared(i,k));

      //apply numerical dampening
      double HrLeft, HrRight, HzBot, HzTop;

      HrLeft = leftFace - Dleft;
      HrRight = rightFace - Dright;
      HzBot = botFace - Dbot;
      HzTop = topFace - Dtop;

      mvRNew(i,k) = mvROld(i,k) - dt*((HrRight-HrLeft)/dr + (HzTop-HzBot)/dz - mvSource(i,k));//update new value
    }//end not halo points
  }//end inside grid
}

__global__ void updateMvZhat(double *mOld, double*mvZNew, double *mvROld, double *mvZOld,
  double *mvSource, int nr, int nz, double dr, double dz, double dt, bool polarity,
  int atomicMass, double rin, double rout, double propellantFlowRate, double a){
  extern __shared__ double U_s[];
  //divide up shared memory
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
  if(polarity){//positive ions
    Mparticle = MASS_POSITIVE_ION;
    Temperature = Tp;
  }
  else{//negative ions
    Mparticle = MASS_NEGATIVE_ION;
    Temperature = Tn;
  }

  //Copy data into shared memory for speedup
    if((i != -1 && i != nr+1)){//not r ghost point
      if(k == -1 || k == nz+1){//z ghost point
        if(k==-1){//thruster inlet
          mvROldShared(0,0)=0;//no radial velocity
          mvZOldShared(0,0)=INLET_MOMENTUM((i+0.5));//Based on propellant flow rate
          mOldShared(0,0)=mvZOldShared(0,0)/a;//Incoming mass
        }
        else{//thruster exit
          //use continuous gradient to approximate
          mvROldShared(0,0)=2*mvROld(i,(k-1))-mvROld(i,(k-2));
          mvZOldShared(0,0) = 2*mvZOld(i,(k-1))-mvZOld(i,(k-2));
          mOldShared(0,0) = 2*massOld(i,(k-1))-massOld(i,(k-2));
        }
      }
      else if(i > -1 && i< nr + 1 && k > -1 && k < nz + 1){//inside grid
        //move to shared memory
        mvROldShared(0,0) = mvROld(i,k);
        mvZOldShared(0,0) = mvZOld(i,k);
        mOldShared(0,0) = massOld(i,k);
      }
    }

  __syncthreads();//Wait for all memory copy to be done

  if(i > -1 && i< nr + 1 && k > -1 && k < nz + 1){//inside gird
    if(threadIdx.x != 0 && threadIdx.x != R_EVALS_PER_BLOCK + 1 && threadIdx.y != 0 && threadIdx.y != Z_EVALS_PER_BLOCK + 1){//not halo points
      //Calculate fluxes
      double botFace, topFace, leftFace, rightFace;
      double Dbot, Dtop, Dleft, Dright;


      if(i==0){//center of thruster
        leftFace = 0;//no streams across the center of thruster
        rightFace = 0.5*((mvZOldShared(0,0)*mvROldShared(0,0)/mOldShared(0,0))+(mvZOldShared(1,0)*mvROldShared(1,0)/mOldShared(1,0)));
        //Numerical dampening
        Dleft = 0;
        Dright = 0.5*0.5*((mvROldShared(1,0)/mOldShared((i+1),k)+a)+(mvROldShared(0,0)/mOldShared(i,k)+a))*(mvZOldShared((i+1),k)-mvZOldShared(i,k));
      }
      else if(i==nr){//top of thruster
        leftFace = 0.5*((mvZOldShared((-1),0)*mvROldShared((-1),0)/mOldShared((-1),0))+(mvZOldShared(0,0)*mvROldShared(0,0)/mOldShared(0,0)));
        rightFace = 0;//no streams out the top of the thruster
        //Numerical dampening
        Dleft = 0.5*0.5*((mvROldShared(0,0)/mOldShared(i,k)+a)+(mvROldShared((-1),0)/mOldShared((i-1),k)+a))*(mvZOldShared(i,k)-mvZOldShared((i-1),k));
        Dright = 0;
      }
      else{//inside thruster
        leftFace = 0.5*((mvZOldShared((-1),0)*mvROldShared((-1),0)/mOldShared((-1),0))+(mvZOldShared(0,0)*mvROldShared(0,0)/mOldShared(0,0)));
        rightFace = 0.5*((mvZOldShared(0,0)*mvROldShared(0,0)/mOldShared(0,0))+(mvZOldShared(1,0)*mvROldShared(1,0)/mOldShared(1,0)));
        //Numerical dampening
        Dleft = 0.5*0.5*((mvROldShared(0,0)/mOldShared(i,k)+a)+(mvROldShared((-1),0)/mOldShared((i-1),k)+a))*(mvZOldShared(i,k)-mvZOldShared((i-1),k));
        Dright = 0.5*0.5*((mvROldShared(1,0)/mOldShared((i+1),k)+a)+(mvROldShared(0,0)/mOldShared(i,k)+a))*(mvZOldShared((i+1),k)-mvZOldShared(i,k));
      }
      botFace = 0.5*((mvZOldShared(0,(-1))*mvZOldShared(0,(-1))/mOldShared(0,(-1))+mOldShared(0,(-1))*kBoltz*Temperature/Mparticle)+(mvZOldShared(0,0)*mvZOldShared(0,0)/mOldShared(0,0)+mOldShared(0,0)*kBoltz*Temperature/Mparticle));
      topFace = 0.5*((mvZOldShared(0,0)*mvZOldShared(0,0)/mOldShared(0,0)+mOldShared(0,0)*kBoltz*Temperature/Mparticle)+(mvZOldShared(0,1)*mvZOldShared(0,1)/mOldShared(0,1)+mOldShared(0,1)*kBoltz*Temperature/Mparticle));
      //Numerical dampening
      Dbot = 0.5*0.5*((mvZOldShared(0,0)/mOldShared(i,k)+a)+(mvZOldShared(0,(-1))/mOldShared(i,(k-1))+a))*(mvZOldShared(i,k)-mvZOldShared(i,(k-1)));
      Dtop = 0.5*0.5*((mvZOldShared(0,1)/mOldShared(i,(k+1))+a)+(mvZOldShared(0,0)/mOldShared(i,k)+a))*(mvZOldShared(i,(k+1))-mvZOldShared(i,k));

      //apply numerical dampening
      double HrLeft, HrRight, HzBot, HzTop;

      HrLeft = leftFace - Dleft;
      HrRight = rightFace - Dright;
      HzBot = botFace - Dbot;
      HzTop = topFace - Dtop;


      mvZNew(i,k) = mvZOld(i,k) - dt*((HrRight-HrLeft)/dr + (HzTop-HzBot)/dz - mvSource(i,k));//update
    }//end not halo points
  }//end inside grid
}

void getMomentum(double *mOldP, double *mvRNewP, double *mvZNewP, double *mvROldP, double *mvZOldP, double *mvRSourceP, double *mvZSourceP,
  double *mOldN, double *mvRNewN, double *mvZNewN, double *mvROldN, double *mvZOldN, double *mvRSourceN, double *mvZSourceN,
  int nr, int nz, double dr, double dz, double dt, int atomicMass, double rin, double rout, double propellantFlowRate, double a,
  dim3 centerGridWHalosBlockDim, dim3 centerGridWHalosThreadDim){

      updateMvRhat<<<centerGridWHalosBlockDim,centerGridWHalosThreadDim,3*(R_EVALS_PER_BLOCK+2)*(Z_EVALS_PER_BLOCK+2)*sizeof(double)>>>
      (mOldP, mvRNewP, mvROldP, mvZOldP, mvRSourceP, nr, nz, dr, dz, dt, true, atomicMass, rin, rout, propellantFlowRate, a);//Update r hat momentum positives

      updateMvRhat<<<centerGridWHalosBlockDim,centerGridWHalosThreadDim,3*(R_EVALS_PER_BLOCK+2)*(Z_EVALS_PER_BLOCK+2)*sizeof(double)>>>
      (mOldN, mvRNewN, mvROldN, mvZOldN, mvRSourceN, nr, nz, dr, dz, dt, false, atomicMass, rin, rout, propellantFlowRate, a);//Update r hat momentum negatives

      updateMvZhat<<<centerGridWHalosBlockDim,centerGridWHalosThreadDim,3*(R_EVALS_PER_BLOCK+2)*(Z_EVALS_PER_BLOCK+2)*sizeof(double)>>>
      (mOldP, mvZNewP, mvROldP, mvZOldP, mvZSourceP, nr, nz, dr, dz, dt, true, atomicMass, rin, rout, propellantFlowRate, a);//Update z hat momentum positives

      updateMvZhat<<<centerGridWHalosBlockDim,centerGridWHalosThreadDim,3*(R_EVALS_PER_BLOCK+2)*(Z_EVALS_PER_BLOCK+2)*sizeof(double)>>>
      (mOldN, mvZNewN, mvROldN, mvZOldN, mvZSourceN, nr, nz, dr, dz, dt, false, atomicMass, rin, rout, propellantFlowRate, a);//Update z hat momentum negatives

  }


#endif
