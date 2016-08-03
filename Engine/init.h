#ifndef _INIT_H_
#define _INIT_H_

__global__ void initCenterQuantities(double *np, double *nn, double *npr,
  double *nnr, double *npz, double *nnz, int nump, int numn, double vpr,
  double vnr, double vpz, double vnz, double dr, double rin, double nr){

  int i = blockIdx.x*R_EVALS_PER_BLOCK + threadIdx.x;
  int k = blockIdx.y*Z_EVALS_PER_BLOCK + threadIdx.y;
  int gridPos = k*(nr+1)+i;
  np[gridPos] = nump*((i+0.5)+rin)*dr;
  nn[gridPos] = numn*((i+0.5)+rin)*dr;
  npr[gridPos] = nump*((i+0.5)+rin)*dr*vpr;
  nnr[gridPos] = numn*((i+0.5)+rin)*dr*vnr;
  npz[gridPos] = nump*((i+0.5)+rin)*dr*vpz;
  nnz[gridPos] = numn*((i+0.5)+rin)*dr*vnz;
}

void getInit(double *np, double *nn, double *npr, double *nnr, double *npz, double *nnz,
  int nump, int numn, double vpr,
  double vnr, double vpz, double vnz, double dr, double rin, double nr,
  dim3 centerGridNoHalosBlockDim, dim3 centerGridNoHalosThreadDim){

  initCenterQuantities<<<centerGridNoHalosBlockDim,centerGridNoHalosThreadDim>>>
    (np,nn,npr,nnr,npz,nnz,nump,numn,vpr,vnr,vpz,vnz,dr,rin,nr);
}

#endif
