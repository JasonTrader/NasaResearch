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

  //Time loop
  double t = startTime;
  while(t < endTime){

    //TODO Calculate electric potential

    //TODO calculate fluxes

    //TODO Update conserved quamtities

    //TODO Update secondary quantities


  }

  //TODO Output results

  //TODO free memory

  return 0;
}
