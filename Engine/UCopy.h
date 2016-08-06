#ifndef _UCOPY_H_
#define _UCOPY_H_

//Copy data from new vectors to old vectors in order to synchronize across
//blocks
void UCopy(double *mOldP,double *mNewP, double *mOldN, double *mNewN,
  double *mvROldP, double *mvRNewP, double *mvROldN, double *mvRNewN,
  double *mvZOldP, double *mvZNewP, double *mvZOldN, double *mvZNewN, size_t centerGridSize){

  cudaMemcpy(mOldP, mNewP, centerGridSize, cudaMemcpyDeviceToDevice);
  cudaMemcpy(mOldN, mNewN, centerGridSize, cudaMemcpyDeviceToDevice);
  cudaMemcpy(mvROldP, mvRNewP, centerGridSize, cudaMemcpyDeviceToDevice);
  cudaMemcpy(mvROldN, mvRNewN, centerGridSize, cudaMemcpyDeviceToDevice);
  cudaMemcpy(mvZOldP, mvZNewP, centerGridSize, cudaMemcpyDeviceToDevice);
  cudaMemcpy(mvZOldN, mvZNewN, centerGridSize, cudaMemcpyDeviceToDevice);
}

#endif
