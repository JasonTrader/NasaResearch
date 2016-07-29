#include <stdlib.h>
#include <stdio.h>

//Header files
#include "LaplacianVoltage.h"

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

//---------------------------------------------------------------------------//

  //TODO calculate initial conserved quantities

  //TODO calculate secondary initial quantities


//---------------------------------------------------------------------------//
//This area is used for thread/block dim calcs

  dim3 cornerGridWHalosThreadDim(R_EVALS_PER_BLOCK+2, Z_EVALS_PER_BLOCK+2);//+ 2 accounts for halo points
  int cornerGridWHalosR = 1 + (nr-1)/R_EVALS_PER_BLOCK;//nr = number internal r grid points
  int cornerGridWHalosZ = 1 + (nz-1)/Z_EVALS_PER_BLOCK;//nz = number internal z grid points
  dim3 cornerGridWHalosBlockDim(cornerGridWHalosR, cornerGridWHalosZ);


  //Time loop
  double t = startTime;
  while(t < endTime){


//---------------------------------------------------------------------------//
//Update Voltage
    bool didConverge = false;
    while(!didConverge){
      cudaMemcpy(volt_d_old, volt_d_new, gridSize, cudaMemcpyDeviceToDevice);
      //Evaluate red blocks
      updateVoltage<<<cornerGridWHalosBlockDim,cornerGridWHalosThreadDim,(zEvalsPerBlock+2)*(rEvalsPerBlock+2)*sizeof(double)>>>(volt_d_old, volt_d_new, true, converge_d, nr, nz, dr, dz, blockz);
      cudaMemcpy(volt_d_old, volt_d_new, gridSize, cudaMemcpyDeviceToDevice);
      //Evaluate black blocks
      updateVoltage<<<cornerGridWHalosBlockDim,cornerGridWHalosThreadDim,(zEvalsPerBlock+2)*(rEvalsPerBlock+2)*sizeof(double)>>>(volt_d_old, volt_d_new, false, converge_d, nr, nz, dr, dz, blockz);
      //copy back converge check
      cudaMemcpy(converge_h, converge_d, convSize, cudaMemcpyDeviceToHost);

      //all converge must be true
      didConverge = converge_h[0];
      for(int i = 1; i< 2*blockr*blockz; i++){
        didConverge = didConverge && converge_h[i];
      }
    }//converged

//---------------------------------------------------------------------------//

    //TODO calculate fluxes

    //TODO Update conserved quamtities

    //TODO Update secondary quantities

    t+=dt;//update time
  }

  //TODO Output results

  //TODO free memory

  return 0;
}
