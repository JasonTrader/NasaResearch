#ifndef _DTFIND_H_
#define _DTFIND_H_

//FIXME speed up move to source sink
__global__ void dtFind(double *posMass, double *negMass, double *momentumRp,
  double *momentumRn, double *momentumZp, double *momentumZn, double *smallestDt,
  double dr, double dz, int nr, int nz, double a, int blockXDim){

  int i = blockIdx.x*R_EVALS_PER_BLOCK + threadIdx.x;
  int k = blockIdx.y*Z_EVALS_PER_BLOCK + threadIdx.y;
  for(int w = 0; w< R_EVALS_PER_BLOCK*Z_EVALS_PER_BLOCK; w++){//all threads
    if(w%R_EVALS_PER_BLOCK==threadIdx.x){//right r
      if(w/R_EVALS_PER_BLOCK==threadIdx.y){//right z
        if(i<nr+1 && k <nz+1){//in grid

          //get largest speed + speed of sound
          double positiveR, negativeR, positiveZ, negativeZ;
          positiveR = fabs(momentumRp[k*(nr+1)+i]/posMass[k*(nr+1)+i]) +fabs(a);
          negativeR = fabs(momentumRn[k*(nr+1)+i]/negMass[k*(nr+1)+i]) +fabs(a);
          positiveZ = fabs(momentumZp[k*(nr+1)+i]/posMass[k*(nr+1)+i]) +fabs(a);
          negativeZ = fabs(momentumZn[k*(nr+1)+i]/negMass[k*(nr+1)+i]) +fabs(a);


          double smallest = (dr < dz ? dr : dz);

          //getDeltaT
          double dtRp, dtRn, dtZp, dtZn;
          dtRp = 0.5*smallest*positiveR;
          dtRn = 0.5*smallest*negativeR;
          dtZp = 0.5*smallest*positiveZ;
          dtZn = 0.5*smallest*negativeZ;


          double min, first, second;
          first = (dtRp < dtRn ? dtRp : dtRn);
          second = (dtZp < dtZn ? dtZp : dtZn);
          min = (first < second ? first : second);
          if(smallestDt[blockIdx.y*blockXDim+blockIdx.x] > min){
            smallestDt[blockIdx.y*blockXDim+blockIdx.x] = min;
          }

        }//end grid
      }//end z
    }//end r
    __syncthreads();
  }//end loop
}

double getDt(double *posMass, double *negMass, double *momentumRp,
  double *momentumRn, double *momentumZp, double *momentumZn,
  double *smallestDt_d, double *smallestDt_h, int blockR, int blockZ,
  double dr, double dz, int nr, int nz, double a,
  dim3 centerGridNoHalosBlockDim, dim3 centerGridNoHalosThreadDim){

    dtFind<<<centerGridNoHalosBlockDim,centerGridNoHalosThreadDim>>>
      (posMass,negMass,momentumRp,momentumRn,momentumZp,momentumZn,smallestDt_d,dr,dz,nr,nz,a,blockR);

    cudaMemcpy(smallestDt_h, smallestDt_d, blockR*blockZ*sizeof(double), cudaMemcpyDeviceToHost);
    //FIXME speed up with sort
    double min = smallestDt_h[0];
    for(int i = 1; i<blockR*blockZ; i++){
      if(min > smallestDt_h[i])
        min = smallestDt_h[i];
    }
    return min;

  }

#endif
