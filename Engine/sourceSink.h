#ifndef _SOURCESINK_H_
#define _SOURCESINK_H_

#include "globals.h"

__device__ double getGamma(){
  return 0;
}

__global__ void updateSourceSink(double *smp, double *smn, double *smvRp, double *smvRn, double *smvZp, double *smvZn,
                                double *mp, double *mn, double *mvRp, double *mvRn, double *mvZp, double *mvZn,
                                double *Er, double *Ez, double dr, int nr){
  int i = blockIdx.x * R_EVALS_PER_BLOCK + threadIdx.x;//x position index
  int k = blockIdx.y * Z_EVALS_PER_BLOCK + threadIdx.y;//y position index
  int gridPos = k*(nr+1)+i;
  int np = mp[gridPos]/((i+0.5)*dr);
  int nn = mn[gridPos]/((i+0.5)*dr);
  double vpr = mvRp[gridPos]/mp[gridPos];
  double vnr = mvRn[gridPos]/mn[gridPos];
  double vpz = mvZp[gridPos]/mp[gridPos];
  double vnz = mvZn[gridPos]/mn[gridPos];

  smp[gridPos] = 0;//Mass positive source
  smn[gridPos] = 0;//Mass negative source
  smvRp[gridPos] = mp[gridPos]*((q*Er[gridPos]/Mp)-(getGamma()*vpr))+np*kBoltz*Tp/Mp;//Radial positive momentum source
  smvRn[gridPos] = mn[gridPos]*((-q*Er[gridPos]/Mn)-(getGamma()*vnr))+nn*kBoltz*Tn/Mn;//Radial negative momentum source
  smvZp[gridPos] = mp[gridPos]*((q*Ez[gridPos]/Mp)-(getGamma()*vpz));//Axial positive momentum source
  smvZn[gridPos] = mn[gridPos]*((-q*Ez[gridPos]/Mn)-(getGamma()*vnz));//Axial negative momentum source
}

void getSourceSink(double *smp, double *smn, double *smvRp, double *smvRn, double *smvZp, double *smvZn,
                    double *mp, double *mn, double *mvRp, double *mvRn, double *mvZp, double *mvZn,
                    double *Er, double *Ez, double dr, int nr, dim3 centerGridNoHalosBlockDim, dim3 centerGridNoHalosThreadDim){
  updateSourceSink<<<centerGridNoHalosBlockDim,centerGridNoHalosThreadDim>>>
                    (smp,smn,smvRp,smvRn,smvZp,smvZn,mp,mn,mvRp,mvRn,mvZp,mvZn,Er,Ez,dr,nr);
}

#endif
