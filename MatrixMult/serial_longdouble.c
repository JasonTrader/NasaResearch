#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void matMul(double *a, double *b, double *c, int arow, int acol, int brow, int bcol){
  int crow = arow;
  int ccol = bcol;
  int i, j;
  for(i = 0; i<crow; i++){
    for(j = 0; j<ccol; j++){
      double sum = 0;
      int k;
      for(k = 0; k< acol; k++){
        sum+=(a[i*acol + k] * b[k*bcol + j]);
      }
      c[i*ccol + j] = sum;
    }
  }
}
int main(){
  double * mat_a, *mat_b, *mat_c;
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
  size_t asiz = arow*acol*sizeof(double);
  size_t bsiz = brow*bcol*sizeof(double);
  size_t csiz = crow*ccol*sizeof(double);
  void *temp_a = malloc(asiz);
  void *temp_b = malloc(bsiz);
  void *temp_c = malloc(csiz);
  mat_a = (double *)temp_a;
  mat_b = (double *)temp_b;
  mat_c = (double *)temp_c;


  //initialize
  //values are row driven
  //meaning mat(row,col) = *(mat + row*colMax + col)

  int i, j;
  for(i = 0; i<arow; i++){
    for(j = 0; j<acol; j++){
      mat_a[i*acol+j] = (double)(rand()%10);
    }
  }

  for(i = 0; i<brow; i++){
    for(j = 0; j<bcol; j++){
      mat_b[i*bcol+j] = (double)(rand()%10);
    }
  }

  //solve a*b
  clock_t begin = clock();
  matMul(mat_a, mat_b, mat_c, arow, acol, brow, bcol);
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  //output results

  FILE *f;
  f = fopen("serRes.txt", "w");
  for(i = 0; i<crow; i++){
    for(j = 0; j<ccol; j++){
//      fprintf(f, "%lf ", mat_c[i*ccol+j]);
      fprintf(f, "%lf ", mat_c[i*ccol+j]);
    }
    fprintf(f, "\n");
  }
  fclose(f);


  //free memory
  free(mat_a);
  free(mat_b);
  free(mat_c);

  printf("\n--------------------------\nExecution took: %lf seconds\n", time_spent);

  return 0;
}
