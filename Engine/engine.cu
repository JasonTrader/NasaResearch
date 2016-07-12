#include <stdlib.h>
#include <stdio.h>

#define PI 3.1415926535897932384
#define mu0 4*PI*1e-7
//Threads per block is capped at 1024 for hardware reasons
//In some cases using a smaller number of threads per block will be more efficient
#define threadsPerBlock 1024
//Max grid points is to defined in order to allocate shared memory for each block
#define MaxGridPoints 6144


int main(){
  int numProblems;
  scanf("%d", &numProblems);

  for(int i =0; i< numProblems; i++){
    //TODO getInitial
  }

  //TODO Copy needed memory down

  //TODO Run Init function

  //TODO Copy back orignal conditions for output maybe?

  //TODO output inital to file maybe?

  //TODO run simulation

  //TODO copy back results

  //TODO Output results

  //TODO free memory
  
  return 0;
}
