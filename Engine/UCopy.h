#ifndef _UCOPY_H_
#define _UCOPY_H_

void UCopy(double *mOldP,double *mNewP, double *mOldN, double *mNewN,
  double *mvROldP, double *mvRNewP, double *mvROldN, double *mvRNewN,
  double *mvZOldP, double *mvZNewP, double *mvZOldN, double *mvZNewN, size_t centerGridSize){

  cudaMemcpy(mNewP, mOldP, centerGridSize, cudaMemcpyDeviceToDevice);
  cudaMemcpy(mNewN, mOldN, centerGridSize, cudaMemcpyDeviceToDevice);
  cudaMemcpy(mvRNewP, mvROldP, centerGridSize, cudaMemcpyDeviceToDevice);
  cudaMemcpy(mvRNewN, mvROldN, centerGridSize, cudaMemcpyDeviceToDevice);
  cudaMemcpy(mvZNewP, mvZOldP, centerGridSize, cudaMemcpyDeviceToDevice);
  cudaMemcpy(mvZNewN, mvZOldN, centerGridSize, cudaMemcpyDeviceToDevice);
}

#endif
