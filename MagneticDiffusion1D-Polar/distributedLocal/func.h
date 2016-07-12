#ifndef FUNC_H
#define FUNC_H

#define PI 3.1415926535897932384
#define mu0 4*PI*1e-7

__global__ void init(double *rod_new, double *rod_old, int *numTimeSteps, bool *ready, double totCurrent, double ldr, double rlength, int numPoints) {
  int i = blockIdx.x;
  if(i==0){
    rod_new[i] = 0;
    rod_old[i] = 0;
  }
  else if (i < (numPoints-1)){
    double init_val = (1-(i*i*ldr*ldr/(3*rlength*rlength)))*3*mu0*totCurrent*i*ldr/(4*PI*rlength*rlength);
    rod_new[i] = init_val;
    rod_old[i] = init_val;
    ready[i] = false;
    numTimeSteps[i] = 0;
  }
  else if(i==(numPoints-1)){
    double init_val = mu0*totCurrent/(2*PI*rlength);
    rod_new[i] = init_val;
    rod_old[i] = init_val;
  }

}

//num points is different than the above mentioned num points
__global__ void run(double *rod_new, double *rod_old, int numPoints, double aug, int maxT, bool *ready, int *numTimeSteps){
  int i = blockIdx.x + 1;
  if (i < numPoints){
    if(i==1){
      while(numTimeSteps[i] <= maxT){
        while(numTimeSteps[i] > numTimeSteps[i+1]){
          //printf("Waiting to move forward at \t%d\t%d\n", i, numTimeSteps[i]);
        }
        rod_old[i] = rod_new[i];
        ready[i] = !ready[i];
        bool check = numTimeSteps[i]%2;
        while((numTimeSteps[i+1]%2==1) != check){
          //printf("Waiting for neighbor to copy\t%d\t%d\n", i, numTimeSteps[i]);
        }
        rod_new[1]+= aug*(2*rod_old[2] - 4*rod_old[1]);
        numTimeSteps[i]++;
      }
    }
    else if(i==(numPoints-1)){
      while(numTimeSteps[i] <= maxT){
        while(numTimeSteps[i] > numTimeSteps[i-1]){
          //printf("Waiting to move forward at \t%d\t%d\n", i, numTimeSteps[i]);
        }
        rod_old[i] = rod_new[i];
        ready[i] = !ready[i];
        bool check = numTimeSteps[i]%2;
        while((numTimeSteps[i-1]%2==1) != check){
          //printf("Waiting for neighbor to copy\t%d\t%d\n", i, numTimeSteps[i]);
        }
        rod_new[i] += aug*((1+(1/(2*i)))*rod_old[i+1] + (-2-(1/(i*i)))*rod_old[i] + (1-(1/(2*i)))*rod_old[i-1]);
        numTimeSteps[i]++;
      }
    }
    else{
      while(numTimeSteps[i] <= maxT){
        while(numTimeSteps[i] > numTimeSteps[i+1] || numTimeSteps[i] > numTimeSteps[i-1]){
          //printf("Waiting to move forward at \t%d\t%d\n", i, numTimeSteps[i]);
        }
        rod_old[i] = rod_new[i];
        ready[i] = !ready[i];
        bool check = numTimeSteps[i]%2;
        while((numTimeSteps[i+1]%2==1) != check || (numTimeSteps[i-1]%2==1) != check){
          //printf("Waiting for neighbor to copy\t%d\t%d\n", i, numTimeSteps[i]);
        }
        rod_new[i] += aug*((1+(1/(2*i)))*rod_old[i+1] + (-2-(1/(i*i)))*rod_old[i] + (1-(1/(2*i)))*rod_old[i-1]);
        numTimeSteps[i]++;
      }//end of while
    }//end of else

    printf("Point #%d at step %d\n", i, numTimeSteps[i]);

  }//end of i<numPoints
}//end of run

#endif
