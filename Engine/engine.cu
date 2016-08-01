#include <stdlib.h>
#include <stdio.h>

//Header files
#include "LaplacianVoltage.h"
#include "MassConservation.h"

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
  double dt = 0.125*smallest*smallest*MU0/eta;//to ensure stability
  //TODO eta?

//---------------------------------------------------------------------------//
// Memory setup

  size_t cornerGridSize = (nr+2)*(nz+2)*sizeof(double);
  dim3 cornerGridWHalosThreadDim(R_EVALS_PER_BLOCK+2, Z_EVALS_PER_BLOCK+2);//+ 2 accounts for halo points
  int cornerGridWHalosBlockR = 1 + (nr-1)/R_EVALS_PER_BLOCK;//nr = number internal r grid points
  int cornerGridWHalosBlockZ = 1 + (nz-1)/Z_EVALS_PER_BLOCK;//nz = number internal z grid points
  dim3 cornerGridWHalosBlockDim(cornerGridWHalosBlockR, cornerGridWHalosBlockZ);

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

    //TODO calculate fluxes

    //NOTE these two steps may be able to be combined

    //TODO Update conserved quantities

    //TODO Update secondary quantities

    t+=dt;//update time
  }

  //TODO Output results

//---------------------------------------------------------------------------//
//Free memory
  cudaFree(voltOld_d);
  cudaFree(voltNew_d);
  free(volt_h);

  return 0;
}
