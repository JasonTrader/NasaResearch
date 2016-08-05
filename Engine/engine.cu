#include <stdlib.h>
#include <stdio.h>

//Header files
#include "globals.h"
#include "params.h"
#include "init.h"
#include "LaplacianVoltage.h"
#include "sourceSink.h"
#include "MassConservation.h"
#include "MomentumConservation.h"
#include "UCopy.h"
#include "dtFind.h"

#define U_d(i,on) U_d+((2*i+on)*(nr+1)*(nz+1))
#define U_h(i) U_h+(i*(nr+1)*(nz+1))
#define S_d(i) S_d+(i*(nr+1)*(nz+1))
#define S_h(i) S_h+(i*(nr+1)*(nz+1))

int main(){

  //-------------------------------------------------------------------------//
  //get input
  int atomicMass, nr, nz;
  double propellantFlowRate,rIn, rOut, lr, lz, startTime, endTime, vMax, vMin,dr,dz,biasF,a;
  getParams(&atomicMass,&propellantFlowRate,&rIn,&rOut,&lr,&lz,&nr,&nz,&startTime,&endTime,&dr,&dz,&vMax,&vMin,&biasF);
  a = (Tp > Tn ? sqrt(2*kBoltz*Tp/MASS_POSITIVE_ION) : sqrt(2*kBoltz*Tn/MASS_POSITIVE_ION));//speed of sound

//---------------------------------------------------------------------------//
// Memory setup

  size_t cornerGridSize = (nr+2)*(nz+2)*sizeof(double);
  dim3 cornerGridWHalosThreadDim(R_EVALS_PER_BLOCK+2, Z_EVALS_PER_BLOCK+2);//+ 2 accounts for halo points
  int cornerGridWHalosBlockR = 1 + (nr-1)/R_EVALS_PER_BLOCK;//nr = number internal r grid points
  int cornerGridWHalosBlockZ = 1 + (nz-1)/Z_EVALS_PER_BLOCK;//nz = number internal z grid points
  dim3 cornerGridWHalosBlockDim(cornerGridWHalosBlockR, cornerGridWHalosBlockZ);

  size_t centerGridSize = (nr+1)*(nz+1)*sizeof(double);
  dim3 centerGridWHalosThreadDim(R_EVALS_PER_BLOCK+2, Z_EVALS_PER_BLOCK+2);//+ 2 accounts for halo points
  int centerGridBlockR = 1 + (nr+1-1)/R_EVALS_PER_BLOCK;//nr+1 = number internal r grid points
  int centerGridBlockZ = 1 + (nz+1-1)/Z_EVALS_PER_BLOCK;//nz + 1 = number internal z grid points
  dim3 centerGridWHalosBlockDim(centerGridBlockR, centerGridBlockZ);
  dim3 centerGridNoHalosThreadDim(R_EVALS_PER_BLOCK, Z_EVALS_PER_BLOCK);
  dim3 centerGridNoHalosBlockDim(centerGridBlockR, centerGridBlockZ);


  //Voltage
  double *voltOld_d, *voltNew_d;//Device voltage grids
  cudaMalloc(&voltOld_d, cornerGridSize);
  cudaMalloc(&voltNew_d, cornerGridSize);
  double *volt_h;//Host voltage grid
  volt_h = (double*)malloc(cornerGridSize);
  bool *converge_d;
  size_t convergeSize = 2*cornerGridWHalosBlockR*cornerGridWHalosBlockZ*sizeof(bool);
  cudaMalloc(&converge_d, convergeSize);
  bool *converge_h;
  converge_h = (bool*)malloc(convergeSize);

  //E Fields
  double *Er_d, *Ez_d;
  cudaMalloc(&Er_d,centerGridSize);
  cudaMalloc(&Ez_d,centerGridSize);
  double *Er_h, *Ez_h;
  Er_h = (double*)malloc(centerGridSize);
  Ez_h = (double*)malloc(centerGridSize);

  //Conserved Quantities
  size_t uSize = 6*centerGridSize;
  double *U_d, *S_d;
  cudaMalloc(&U_d,2*uSize);
  cudaMalloc(&S_d,uSize);
  double *U_h, *S_h;
  U_h = (double*)malloc(centerGridSize);
  S_h = (double*)malloc(centerGridSize);

  //dt
  double *dtVec_d;
  cudaMalloc(&dtVec_d,centerGridBlockR*centerGridBlockZ*sizeof(double));
  double * dtVec_h;
  dtVec_h = (double*)malloc(centerGridBlockR*centerGridBlockZ*sizeof(double));


//---------------------------------------------------------------------------//
  double np = 1e14;//#/m^3
  double nn = 1e14;//#/m^3
  double Vpr = 0;//Initial velocity
  double Vnr = 0;//Initial velocity
  double Vpz = 0;//Initial velocity
  double Vnz = 0;//Initial velocity

  getInit(U_d(massP,o),U_d(massN,o),U_d(momentumPR,o),U_d(momentumNR,o),U_d(momentumPZ,o),U_d(momentumNZ,o),
    np,nn,Vpr,Vnr,Vpz,Vnz,dr,rIn,nr,centerGridNoHalosBlockDim,centerGridNoHalosThreadDim);
  //Time loop
  double t = startTime;
  while(t < endTime){
    double dt = getDt(U_d(massP,o), U_d(massN,o), U_d(momentumPR,o),
      U_d(momentumNR,o), U_d(momentumPZ,o), U_d(momentumNZ,o),
      dtVec_d, dtVec_h, centerGridBlockR, centerGridBlockZ,
      dr,dz,nr,nz,a,centerGridNoHalosBlockDim,centerGridNoHalosThreadDim);


//---------------------------------------------------------------------------//
    //Update Voltage
    getNewVoltage(cornerGridSize,convergeSize,voltOld_d,voltNew_d,U_d(massP,o),U_d(massN,o),volt_h,cornerGridWHalosBlockDim,
      cornerGridWHalosThreadDim,converge_d,converge_h,nr,nz,dr,dz,rIn,vMax,vMin,biasF,t,
      cornerGridWHalosBlockR,cornerGridWHalosBlockZ);

//---------------------------------------------------------------------------//
    getSourceSink(S_d(massP),S_d(massN),S_d(momentumPR),S_d(momentumNR),S_d(momentumPZ),S_d(momentumNZ),
      U_d(massP,o),U_d(massN,o),U_d(momentumPR,o),U_d(momentumNR,o),U_d(momentumPZ,o),U_d(momentumNZ,o),
      Er_d, Ez_d,dr,nr,atomicMass,centerGridNoHalosBlockDim,centerGridNoHalosThreadDim);

    getMass(U_d(massP,o),U_d(massP,n),U_d(momentumPR,o),U_d(momentumPZ,o),S_d(massP),
      U_d(massN,o),U_d(massN,n),U_d(momentumNR,o),U_d(momentumNZ,o),S_d(massN),
      nr,nz,dr,dz,dt,atomicMass,propellantFlowRate,rIn,rOut,a,centerGridWHalosBlockDim,centerGridWHalosThreadDim);

    getMomentum(U_d(massP,o), U_d(momentumPR,n), U_d(momentumPZ,n), U_d(momentumPR,o), U_d(momentumPZ,o), S_d(momentumPR), S_d(momentumPZ),
      U_d(massN,o), U_d(momentumNR,n), U_d(momentumNZ,n), U_d(momentumNR,o), U_d(momentumNZ,o), S_d(momentumNR), S_d(momentumNZ),
      nr,nz,dr,dz,dt,atomicMass,rIn,rOut,propellantFlowRate,a,centerGridWHalosBlockDim,centerGridWHalosThreadDim);

    UCopy(U_d(massP,o),U_d(massP,n), U_d(massN,o), U_d(massN,n),
      U_d(momentumPR,o), U_d(momentumPR,n), U_d(momentumNR,o), U_d(momentumNR,n),
      U_d(momentumPZ,o), U_d(momentumPZ,n), U_d(momentumNZ,o), U_d(momentumNZ,n), centerGridSize);


    t+=dt;//update time
    //TODO Output results
  }


//---------------------------------------------------------------------------//
//Free memory
  cudaFree(voltOld_d);
  cudaFree(voltNew_d);
  free(volt_h);
  cudaFree(Er_d);
  cudaFree(Ez_d);
  free(Er_h);
  free(Ez_h);
  cudaFree(U_d);
  cudaFree(S_d);
  free(U_h);
  free(S_h);

  return 0;
}
