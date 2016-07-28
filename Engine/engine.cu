#include <stdlib.h>
#include <stdio.h>

#define PI 3.1415926535897932384
#define mu0 4*PI*1e-7

int main(){

  //-------------------------------------------------------------------------//
  //get input

  //Propellant (AMU) (currently not used)
  scanf("%*s %*s %*s %*d");

  //Mass flow rate (kg/s) (currently not used)
  scanf("%*s %*s %*s %*lf");

  //inner r (m)
  double rIn;
  scanf("%*s %*s %*s %lf", &rin);

  //outer r (m)
  double rOut;
  scanf("%*s %*s %*s %lf", &rout);

  //Total rlength
  double lr = rout - rin;

  //z length (m)
  double lz;
  scanf("%*s %*s %*s %lf", &lz);

  //number of points in r direction
  int nr;
  scanf("%*s %*s %*s %d", &nr);

  //number of points in z direction
  int nz;
  scanf("%*s %*s %*s %d", &nz);

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
  while(t < totTime){

    //TODO Calculate electric potential

    //TODO calculate fluxes

    //TODO Update conserved quamtities

    //TODO Update secondary quantities


  }

  //TODO Output results

  //TODO free memory

  return 0;
}
