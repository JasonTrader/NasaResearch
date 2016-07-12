#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
using namespace std;

long matMul(long long int *a, long long int *b, long long int *c, int arow, int acol, int brow, int bcol){
  clock_t begin = clock();
  int crow = arow;
  int ccol = bcol;
  int i, j;
  for(i = 0; i<crow; i++){
    for(j = 0; j<ccol; j++){
      int sum = 0;
      int k;
      for(k = 0; k< acol; k++){
        sum+=(a[i*acol + k] * b[k*bcol + j]);
      }
      c[i*ccol + j] = sum;
    }
  }
  clock_t end = clock();
  long time_spent = (long)(end - begin);
  return time_spent;
}

int main(){
  long long int * mat_a, *mat_b, *mat_c;
  int arow, acol, brow, bcol, crow, ccol;
  printf("Dimensions of a: ");
  scanf("%d %d", &arow, &acol);
  printf("Dimensions of b: ");
  scanf("%d %d", &brow, &bcol);

  //check dimensions
  if(acol != brow){
    printf("These matricies may not be multplied together");
    return 1;
  }

  //c's dimensions are a result of multiplication
  crow = arow;
  ccol = bcol;

  //allocate memory
  mat_a = new long long int[arow*acol];
  mat_b = new long long int[brow*bcol];
  mat_c = new long long int[crow*ccol];


  //initialize
  //values are row driven
  //meaning mat(row,col) = *(mat + row*colMax + col)
  int i, j;
  for(i = 0; i<arow; i++){
    for(j = 0; j<acol; j++){
      mat_a[i*acol+j] = i*acol+j;
    }
  }

  for(i = 0; i<brow; i++){
    for(j = 0; j<bcol; j++){
      mat_b[i*bcol+j] = i*bcol+j;
    }
  }

  //solve a*b
  long time_spent = matMul(mat_a, mat_b, mat_c, arow, acol, brow, bcol);

  //output results

  ofstream f;
  f.open("serRes.txt");
  for(i = 0; i<crow; i++){
    for(j = 0; j<ccol; j++){
      f << mat_c[i*ccol+j] << " ";
    }
    f << "\n";
  }
  f.close();


  //free memory
  delete(mat_a);
  delete(mat_b);
  delete(mat_c);

  printf("\n--------------------------\nExecution took: %d clock cycles\n", time_spent);

  return 0;
}
