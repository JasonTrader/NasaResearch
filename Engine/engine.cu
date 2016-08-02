#include <stdlib.h>
#include <stdio.h>

//Header files
#include "globals.h"
#include "LaplacianVoltage.h"
#include "MassConservation.h"
#include "MomentumConservation.h"
#include "UCopy.h"

#define U_d(i,on) U_d+((2*i+on)*(nr+1)*(nz+1))
#define U_h(i) U_h+(i*(nr+1)*(nz+1))
#define S_d(i) S_d+(i*(nr+1)*(nz+1))
#define S_h(i) S_h+(i*(nr+1)*(nz+1))

int main(){

  //-------------------------------------------------------------------------//
  //get input

  //Propellant (AMU) (currently not used)
  scanf("%*s %*s %*s %*d");

  //Mass flow rate (kg/s) (currently not used)
  scanf("%*s %*s %*s %*lf");

  //inner r (m)
  double rIn;
  scanf("%*s %*s %*s %lf", &rIn);

  //outer r (m)
  double rOut;
  scanf("%*s %*s %*s %lf", &rOut);

  //Total rlength
  double lr = rIn - rOut;

  //z length (m)
  double lz;
  scanf("%*s %*s %*s %lf", &lz);

  //number of points in r direction
  int nr;
  scanf("%*s %*s %d", &nr);

  //number of points in z direction
  int nz;
  scanf("%*s %*s %d", &nz);

  //Start time (s)
  double startTime;
  scanf("%*s %*s %*s %lf", &startTime);

  //End time (s)
  double endTime;
  scanf("%*s %*s %*s %lf", &endTime);

//----------------------------------------------------------------------------//
  //variable setup
  double dr = lr/(nr+1);
  double dz = lz/(nz+1);
  double smallest = dr;
  if(dz < dr)
    smallest = dz;
  double dt = 0.125*smallest*smallest*MU0;///eta;//to ensure stability
  //QUESTION eta?

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


//---------------------------------------------------------------------------//

  //TODO calculate initial conserved quantities

  //TODO calculate secondary initial quantities


  //Time loop
  double t = startTime;
  while(t < endTime){


//---------------------------------------------------------------------------//
    //Update Voltage
    getNewVoltage(cornerGridSize,convergeSize,voltOld_d,voltNew_d,volt_h,cornerGridWHalosBlockDim,
      cornerGridWHalosThreadDim,converge_d,converge_h,nr,nz,dr,dz,
      cornerGridWHalosBlockR,cornerGridWHalosBlockZ);

//---------------------------------------------------------------------------//

    getMass(U_d(massP,o),U_d(massP,n),U_d(momentumPR,o),U_d(momentumPZ,o),S_d(massP),
      U_d(massN,o),U_d(massN,n),U_d(momentumNR,o),U_d(momentumNZ,o),S_d(massN),
      nr,nz,dr,dz,dt,centerGridWHalosBlockDim,centerGridWHalosThreadDim,centerGridSize);

    getMomentum(U_d(massP,o), U_d(momentumPR,n), U_d(momentumPZ,n), U_d(momentumPR,o), U_d(momentumPZ,o), S_d(momentumPR), S_d(momentumPZ),
      U_d(massN,o), U_d(momentumNR,n), U_d(momentumNZ,n), U_d(momentumNR,o), U_d(momentumNZ,o), S_d(momentumNR), S_d(momentumNZ),
      nr,nz,dr,dz,dt,centerGridWHalosBlockDim,centerGridWHalosThreadDim,centerGridSize);

    UCopy(U_d(massP,o),U_d(massP,n), U_d(massN,o), U_d(massN,n),
      U_d(momentumPR,o), U_d(momentumPR,n), U_d(momentumNR,o), U_d(momentumNR,n),
      U_d(momentumPZ,o), U_d(momentumPZ,n), U_d(momentumNZ,o), U_d(momentumNZ,n), centerGridSize);

    //TODO Source/Sink

    t+=dt;//update time
  }

  //TODO Output results

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
