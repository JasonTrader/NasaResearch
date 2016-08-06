#ifndef _INIT_H_
#define _INIT_H_

__global__ void initCenterQuantities(double *rnp, double *rnn, double *rnpvpr,
  double *rnnvnr, double *rnpvpz, double *rnnvnz, double nump, double numn, double vpr,
  double vnr, double vpz, double vnz, double dr, double rin, double nr){

  int i = blockIdx.x*R_EVALS_PER_BLOCK + threadIdx.x;//i represents r
  int k = blockIdx.y*Z_EVALS_PER_BLOCK + threadIdx.y;//k represents z
  int gridPos = k*(nr+1)+i;
  rnp[gridPos] = nump*(r((i+0.5)));//rnp
  rnn[gridPos] = numn*(r((i+0.5)));//rnn
  rnpvpr[gridPos] = nump*(r((i+0.5)))*vpr;//rnpvr
  rnnvnr[gridPos] = numn*(r((i+0.5)))*vnr;//rnnvr
  rnpvpz[gridPos] = nump*(r((i+0.5)))*vpz;//rnpvz
  rnnvnz[gridPos] = numn*(r((i+0.5)))*vnz;//rnnvz

  if(i==0 && k==1){
    printf("%e\n", rnp[gridPos]);
  }
//  __syncthreads();
}


void getInit(double *rnp, double *rnn, double *rnpvpr, double *rnnvnr, double *rnpvpz, double *rnnvnz,
  double nump, double numn, double vpr,
  double vnr, double vpz, double vnz, double dr, double rin, double nr,
  dim3 centerGridNoHalosBlockDim, dim3 centerGridNoHalosThreadDim){

  initCenterQuantities<<<centerGridNoHalosBlockDim,centerGridNoHalosThreadDim>>>
    (rnp,rnn,rnpvpr,rnnvnr,rnpvpz,rnnvnz,nump,numn,vpr,vnr,vpz,vnz,dr,rin,nr);
}

#endif
