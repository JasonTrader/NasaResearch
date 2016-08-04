#ifndef _INIT_H_
#define _INIT_H_

__global__ void initCenterQuantities(double *np, double *nn, double *npr,
  double *nnr, double *npz, double *nnz, int nump, int numn, double vpr,
  double vnr, double vpz, double vnz, double dr, double rin, double nr){

  int i = blockIdx.x*R_EVALS_PER_BLOCK + threadIdx.x;//i represents r
  int k = blockIdx.y*Z_EVALS_PER_BLOCK + threadIdx.y;//k represents z
  int gridPos = k*(nr+1)+i;
  np[gridPos] = nump*((i+0.5)+rin)*dr;//rnp
  nn[gridPos] = numn*((i+0.5)+rin)*dr;//rnn
  npr[gridPos] = nump*((i+0.5)+rin)*dr*vpr;//rnpvr
  nnr[gridPos] = numn*((i+0.5)+rin)*dr*vnr;//rnnvr
  npz[gridPos] = nump*((i+0.5)+rin)*dr*vpz;//rnpvz
  nnz[gridPos] = numn*((i+0.5)+rin)*dr*vnz;//rnnvz
}

void getInit(double *np, double *nn, double *npr, double *nnr, double *npz, double *nnz,
  int nump, int numn, double vpr,
  double vnr, double vpz, double vnz, double dr, double rin, double nr,
  dim3 centerGridNoHalosBlockDim, dim3 centerGridNoHalosThreadDim){

  initCenterQuantities<<<centerGridNoHalosBlockDim,centerGridNoHalosThreadDim>>>
    (np,nn,npr,nnr,npz,nnz,nump,numn,vpr,vnr,vpz,vnz,dr,rin,nr);
}

#endif
