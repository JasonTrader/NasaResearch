#ifndef _MASSCONSERVATION_H_
#define _MASSCONSERVATION_H_
//TODO fix r ratio
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
#define mShared(xOffset,yOffset) m_s[(threadIdx.y+yOffset+1)*(R_EVALS_PER_BLOCK+2) + (threadIdx.x+xOffset+1)]

__global__ void updateMass(double *mOld, double*mNew, double *mvR, double *mvZ,
  double *massSource, int nr, int nz, double dr, double dz, double dt,
  int atomicMass, double propellantFlowRate, double rin, double rout, double a){

  extern __shared__ double mv_s[];
  //divide mv_s in halh because CUDA only allows 1 shared vector
  double *mvR_s, *mvZ_s, *m_s;
  mvR_s = mv_s;
  mvZ_s = mv_s + ((R_EVALS_PER_BLOCK+2)*(Z_EVALS_PER_BLOCK+2));
  m_s = mv_s + 2*((R_EVALS_PER_BLOCK+2)*(Z_EVALS_PER_BLOCK+2));
  //in this case when i or k = 0 it is not a boundary position rather
  //an internal volume next to the boundary
  //to get a boundary volume we subtract 1 to enable the creation of a -1
  //"ghost point"
  int i = blockIdx.x*R_EVALS_PER_BLOCK + threadIdx.x - 1;
  int k = blockIdx.y*Z_EVALS_PER_BLOCK + threadIdx.y - 1;

//Copy data into shared memory for speedup
//Note points where i is -1 or nr+1 they are not needed because they are out
//of the scope of the simulation.
  if((i != -1 && i != nr+1)){// not r ghost point
    if(k == -1 || k == nz+1){//z ghost point
      if(k==-1){//Inlet of thruster
        //Can be changed later to accomodate difference of 2e mass
        mvRShared(0,0)=0;//NO r direction velocity at inlet
        mvZShared(0,0)=INLET_MOMENTUM((i+0.5));//Dependent on flow rate
        mShared(0,0)=mvZShared(0,0)/a;
      }
      else{//Thruster exit
        //Use continuous gradient to approximate
        mvRShared(0,0)=2*mvR(i,(k-1))-mvR(i,(k-2));
        mvZShared(0,0)=2*mvZ(i,(k-1))-mvZ(i,(k-2));
        mShared(0,0)=2*mOld(i,(k-1))-mOld(i,(k-2));
      }
    }
    else if(i > -1 && i< nr + 1 && k > -1 && k < nz + 1){//inside grid
      //Move value to shared memory
      mvRShared(0,0) = mvR(i,k);
      mvZShared(0,0) = mvZ(i,k);
      mShared(0,0) = mOld(i,k);
    }
  }

  __syncthreads();//Wait for all memory copy to be done

  if(i > -1 && i< nr + 1 && k > -1 && k < nz + 1){//inside gird
    if(threadIdx.x != 0 && threadIdx.x != R_EVALS_PER_BLOCK + 1 && threadIdx.y != 0 && threadIdx.y != Z_EVALS_PER_BLOCK + 1){//not halo points
      //Calculate fluxes
      double mvZBot, mvZTop, mvRLeft, mvRRight;
      double Dbot, Dtop, Dleft, Dright;

      if(i==0){//Center of rocket
        mvRLeft = 0;//No streams crossing center because of azimuthal consistency
        mvRRight = (r((i+1)))*0.5*(mvRShared(0,0)/(r((i+0.5))) + mvRShared(1,0)/(r((i+1.5))));
        //Numerical dampening
        Dleft = 0;
        Dright = 0.5*0.5*((mvRShared(1,0)/mShared(1,0)+a)+(mvRShared(0,0)/mShared(0,0)+a))*(r((i+1)))*(mShared(1,0)/(r((i+1.5)))-mShared(0,0)/(r((i+0.5))));
      }
      else if(i==nr){//Top of rocket
        mvRLeft = (r(i))*0.5*(mvRShared((-1),0)/(r((i-0.5))) + mvRShared(0,0)/(r((i+0.5))));
        mvRRight = 0;//No streams exiting rocket through the top
        //Numerical dampening
        Dleft = 0.5*0.5*((mvRShared(0,0)/mShared(0,0)+a)+(mvRShared((-1),0)/mShared((-1),0)+a))*(r(i))*(mShared(0,0)/(r((i+0.5)))-mShared((-1),0)/(r((i-0.5))));
        Dright = 0;
      }
      else{//Inside rocket
        mvRLeft = (r(i))*0.5*(mvRShared((-1),0)/(r((i-0.5))) + mvRShared(0,0)/(r((i+0.5))));
        mvRRight = (r((i+1)))*0.5*(mvRShared(0,0)/(r((i+0.5))) + mvRShared(1,0)/(r((i+1.5))));
        //Numerical dampening
        Dleft = 0.5*0.5*((mvRShared(0,0)/mShared(0,0)+a)+(mvRShared((-1),0)/mShared((-1),0)+a))*(r(i))*(mShared(0,0)/(r((i+0.5)))-mShared((-1),0)/(r((i-0.5))));
        Dright = 0.5*0.5*((mvRShared(1,0)/mShared(1,0)+a)+(mvRShared(0,0)/mShared(0,0)+a))*(r((i+1)))*(mShared(1,0)/(r((i+1.5)))-mShared(0,0)/(r((i+0.5))));
      }

      //NOTE r is the same when considering change in z so /r anf *rface is not needed
      //TODO add for not square
      mvZBot = 0.5*(mvZShared(0,(-1)) + mvZShared(0,0));
      mvZTop = 0.5*(mvZShared(0,0) + mvZShared(0,1));
      //Numerical dampening
      Dbot = 0.5*0.5*((mvZShared(0,0)/mShared(0,0)+a)+(mvZShared(0,(-1))/mShared(0,(-1))+a))*(mShared(0,0)-mShared(0,(-1)));
      Dtop = 0.5*0.5*((mvZShared(0,1)/mShared(0,1)+a)+(mvZShared(0,0)/mShared(0,0)+a))*(mShared(0,1)-mShared(0,0));

      //apply numerical dampening
      double HrLeft, HrRight, HzBot, HzTop;

      HrLeft = mvRLeft - Dleft;
      HrRight = mvRRight - Dright;
      HzBot = mvZBot - Dbot;
      HzTop = mvZTop - Dtop;

      if(i==0 && k==0){
        printf("D=%e, delta(r*np)=%e, r*np_Top=%e, r*np_Bot=%e n*np_Old=%e\n",Dtop, (mShared(0,1)-mShared(0,0)), mShared(0,1), mShared(0,0), mOld(0,1));
        //printf("%e %e\n", HzTop, HzBot);
      }

      mNew(i,k) = mShared(0,0) - dt*((HrRight-HrLeft)/dr + (HzTop-HzBot)/dz - massSource(i,k));//Calculates value for next time step
      //if(i==0 && k==0)
        //printf("%e %e\n", mNew(i,k), mOld(i,k));
    }//end not halo points
  }//end inside grid
}

//Calls for both positives and negatives
void getMass(double *mOldP, double *mNewP, double *mvRP, double *mvZP, double *massSourceP,
  double *mOldN, double * mNewN, double *mvRN, double *mvZN, double *massSourceN,
  int nr, int nz, double dr, double dz, double dt, int atomicMass,
  double propellantFlowRate, double rin, double rout, double a,
  dim3 centerGridWHalosBlockDim, dim3 centerGridWHalosThreadDim){

      updateMass<<<centerGridWHalosBlockDim,centerGridWHalosThreadDim,3*(R_EVALS_PER_BLOCK+2)*(Z_EVALS_PER_BLOCK+2)*sizeof(double)>>>
      (mOldP, mNewP, mvRP, mvZP, massSourceP, nr, nz, dr, dz, dt, atomicMass, propellantFlowRate, rin, rout, a);//Update mass positives

      /*updateMass<<<centerGridWHalosBlockDim,centerGridWHalosThreadDim,3*(R_EVALS_PER_BLOCK+2)*(Z_EVALS_PER_BLOCK+2)*sizeof(double)>>>
      (mOldN, mNewN, mvRN, mvZN, massSourceN, nr, nz, dr, dz, dt, atomicMass, propellantFlowRate, rin, rout, a);//Update mass negatives*/
  }


#endif
