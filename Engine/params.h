#ifndef _PARAMS_H_
#define _PARAMS_H_

void getParams(int *atomicMass, double *rIn, double *rOut, double *lr, double *lz,
   int *nr, int *nz, double *startTime, double *endTime, double *dr, double *dz, double *vMax, double *vMin, double *biasF){

  //Propellant (AMU)
  scanf("%*s %*s %*s %d", atomicMass);

  //Mass flow rate (kg/s) (currently not used)
  scanf("%*s %*s %*s %*lf");

  //inner r (m)
  scanf("%*s %*s %*s %lf", rIn);

  //outer r (m)
  scanf("%*s %*s %*s %lf", rOut);

  //Total rlength
  *lr = *rOut - *rIn;

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

  //Vmax (V)
  scanf("%*s %*s %*s %lf", vMax);

  //Vmin (V)
  scanf("%*s %*s %*s %lf", vMin);

  //Bias Frequency (Hz)
  scanf("%*s %*s %*s %lf", biasF);


  *dr = *lr/(*nr+1);
  *dz = *lz/(*nz+1);
}

#endif
