#ifndef _INIT_H_
#define _INIT_H_

void getInit(int *atomicMass, double *rIn, double *rOut, double *lr, double *lz, int *nr, int *nz, double *startTime, double *endTime){
  //TODO test this
  //Propellant (AMU)
  scanf("%*s %*s %*s %d", atomicMass);

  //Mass flow rate (kg/s) (currently not used)
  scanf("%*s %*s %*s %*lf");

  //inner r (m)
  scanf("%*s %*s %*s %lf", rIn);

  //outer r (m)
  scanf("%*s %*s %*s %lf", rOut);

  //Total rlength
  *lr = *rIn - *rOut;

  //z length (m)
  scanf("%*s %*s %*s %lf", lz);

  //number of points in r direction
  scanf("%*s %*s %d", nr);

  //number of points in z direction
  scanf("%*s %*s %d", nz);

  //Start time (s)
  scanf("%*s %*s %*s %lf", startTime);

  //End time (s)
  scanf("%*s %*s %*s %lf", endTime);
}

#endif
