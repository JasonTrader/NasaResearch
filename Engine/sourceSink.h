#ifndef _SOURCESINK_H_
#define _SOURCESINK_H_

#include "globals.h"


__device__ double getNu(){
  return 3e5;
}

__global__ void updateSourceSink(double *smp, double *smn, double *smvRp, double *smvRn, double *smvZp, double *smvZn,
                                double *mp, double *mn, double *mvRp, double *mvRn, double *mvZp, double *mvZn,
                                double *Er, double *Ez, double dr, int nr, int atomicMass){
  int i = blockIdx.x * R_EVALS_PER_BLOCK + threadIdx.x;//x position index
  int k = blockIdx.y * Z_EVALS_PER_BLOCK + threadIdx.y;//y position index
  int gridPos = k*(nr+1)+i;

  int np = mp[gridPos]/((i+0.5)*dr);//rnp
  int nn = mn[gridPos]/((i+0.5)*dr);//rnn
  double vpr = mvRp[gridPos]/mp[gridPos];//rnpvr/rnp
  double vnr = mvRn[gridPos]/mn[gridPos];//rnnvr/rnn
  double vpz = mvZp[gridPos]/mp[gridPos];//rnpvz/rnp
  double vnz = mvZn[gridPos]/mn[gridPos];//rnnvz/rnn

  smp[gridPos] = 0;//Mass positive source
  smn[gridPos] = 0;//Mass negative source
  smvRp[gridPos] = mp[gridPos]*((q*Er[gridPos]/MASS_POSITIVE_ION)-(getNu()*vpr))+np*kBoltz*Tp/MASS_POSITIVE_ION;//Radial positive momentum source
  smvRn[gridPos] = mn[gridPos]*((-q*Er[gridPos]/MASS_NEGATIVE_ION)-(getNu()*vnr))+nn*kBoltz*Tn/MASS_NEGATIVE_ION;//Radial negative momentum source
  smvZp[gridPos] = mp[gridPos]*((q*Ez[gridPos]/MASS_POSITIVE_ION)-(getNu()*vpz));//Axial positive momentum source
  smvZn[gridPos] = mn[gridPos]*((-q*Ez[gridPos]/MASS_NEGATIVE_ION)-(getNu()*vnz));//Axial negative momentum source
}

//Update next steps source sink
void getSourceSink(double *smp, double *smn, double *smvRp, double *smvRn, double *smvZp, double *smvZn,
                    double *mp, double *mn, double *mvRp, double *mvRn, double *mvZp, double *mvZn,
                    double *Er, double *Ez, double dr, int nr, int atomicMass, dim3 centerGridNoHalosBlockDim, dim3 centerGridNoHalosThreadDim){
  updateSourceSink<<<centerGridNoHalosBlockDim,centerGridNoHalosThreadDim>>>
                    (smp,smn,smvRp,smvRn,smvZp,smvZn,mp,mn,mvRp,mvRn,mvZp,mvZn,Er,Ez,dr,nr,atomicMass);
}

#endif
