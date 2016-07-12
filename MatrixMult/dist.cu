#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Thread block size
#define BLOCK_SIZE 16

// Matrices are stored in row-major order:
// M(row, col) = *(M.elements + row * M.stride + col)
typedef struct {
    int width;
    int height;
    int stride;
    long long int* elements;
} Matrix;

// Get a matrix element
__device__ long long int GetElement(const Matrix A, int row, int col)
{
    return A.elements[row * A.stride + col];
}

// Set a matrix element
__device__ void SetElement(Matrix A, int row, int col,
                           long long int value)
{
    A.elements[row * A.stride + col] = value;
}

// Get the BLOCK_SIZExBLOCK_SIZE sub-matrix Asub of A that is
// located col sub-matrices to the right and row sub-matrices down
// from the upper-left corner of A
 __device__ Matrix GetSubMatrix(Matrix A, int blockrow, int blockcol)
{
    Matrix Asub;
    Asub.width    = BLOCK_SIZE;
    Asub.height   = BLOCK_SIZE;
    Asub.stride   = A.stride;
    Asub.elements = &A.elements[A.stride * BLOCK_SIZE * blockrow
                                         + BLOCK_SIZE * blockcol];
    return Asub;
}

// Forward declaration of the matrix multiplication kernel
__global__ void MatMulKernel(const Matrix, const Matrix, Matrix);

// Matrix multiplication - Host code
// Matrix dimensions are assumed to be multiples of BLOCK_SIZE
void MatMul(Matrix A, Matrix B, Matrix C)
{
    // Load A and B to device memory
    Matrix d_A;
    d_A.width = d_A.stride = A.width; d_A.height = A.height;
    size_t size = A.width * A.height * sizeof(long long int);
    cudaMalloc(&d_A.elements, size);
    cudaMemcpy(d_A.elements, A.elements, size,
               cudaMemcpyHostToDevice);
    Matrix d_B;
    d_B.width = d_B.stride = B.width; d_B.height = B.height;
    size = B.width * B.height * sizeof(long long int);
    cudaMalloc(&d_B.elements, size);
    cudaMemcpy(d_B.elements, B.elements, size,
    cudaMemcpyHostToDevice);

    // Allocate C in device memory
    Matrix d_C;
    d_C.width = d_C.stride = C.width; d_C.height = C.height;
    size = C.width * C.height * sizeof(long long int);
    cudaMalloc(&d_C.elements, size);

    // Invoke kernel
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    dim3 dimGrid(B.width / dimBlock.x, A.height / dimBlock.y);
    MatMulKernel<<<dimGrid, dimBlock>>>(d_A, d_B, d_C);

    // Read C from device memory
    cudaMemcpy(C.elements, d_C.elements, size,
               cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(d_A.elements);
    cudaFree(d_B.elements);
    cudaFree(d_C.elements);
}

// Matrix multiplication kernel called by MatMul()
 __global__ void MatMulKernel(Matrix A, Matrix B, Matrix C)
{
    // Block row and column
    int blockRow = blockIdx.y;
    int blockCol = blockIdx.x;

    // Each thread block computes one sub-matrix Csub of C
    Matrix Csub = GetSubMatrix(C, blockRow, blockCol);

    // Each thread computes one element of Csub
    // by accumulating results into Cvalue
    long long int Cvalue = 0;

    // Thread row and column within Csub
    int row = threadIdx.y;
    int col = threadIdx.x;

    // Loop over all the sub-matrices of A and B that are
    // required to compute Csub
    // Multiply each pair of sub-matrices together
    // and accumulate the results
    for (int m = 0; m < (A.width / BLOCK_SIZE); ++m) {

        // Get sub-matrix Asub of A
        Matrix Asub = GetSubMatrix(A, blockRow, m);

        // Get sub-matrix Bsub of B
        Matrix Bsub = GetSubMatrix(B, m, blockCol);

        // Shared memory used to store Asub and Bsub respectively
        __shared__ long long int As[BLOCK_SIZE][BLOCK_SIZE];
        __shared__ long long int Bs[BLOCK_SIZE][BLOCK_SIZE];

        // Load Asub and Bsub from device memory to shared memory
        // Each thread loads one element of each sub-matrix
        As[row][col] = GetElement(Asub, row, col);
        Bs[row][col] = GetElement(Bsub, row, col);

        // Synchronize to make sure the sub-matrices are loaded
        // before starting the computation
        __syncthreads();
        // Multiply Asub and Bsub together
        for (int e = 0; e < BLOCK_SIZE; ++e)
            Cvalue += As[row][e] * Bs[e][col];

        // Synchronize to make sure that the preceding
        // computation is done before loading two new
        // sub-matrices of A and B in the next iteration
        __syncthreads();
    }

    // Write Csub to device memory
    // Each thread writes one element
    SetElement(Csub, row, col, Cvalue);
}

int main(){
  Matrix mat_a, mat_b, mat_c;
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
  size_t asiz = arow*acol*sizeof(long long int);
  size_t bsiz = brow*bcol*sizeof(long long int);
  size_t csiz = crow*ccol*sizeof(long long int);
  void *temp_a = malloc(asiz);
  void *temp_b = malloc(bsiz);
  void *temp_c = malloc(csiz);
  mat_a.elements = (long long int *)temp_a;
  mat_b.elements = (long long int *)temp_b;
  mat_c.elements = (long long int *)temp_c;


  //initialize
  //values are row driven
  //meaning mat(row,col) = *(mat + row*colMax + col)
  long long int i, j;
  for(i = 0; i<arow; i++){
    for(j = 0; j<acol; j++){
      mat_a.elements[i*acol+j] = i*acol+j;
    }
  }

  for(i = 0; i<brow; i++){
    for(j = 0; j<bcol; j++){
      mat_b.elements[i*bcol+j] = i*bcol+j;
    }
  }
  mat_a.height = arow;
  mat_a.width = acol;
  mat_a.stride = acol;
  mat_b.height = brow;
  mat_b.width = bcol;
  mat_b.stride = bcol;
  mat_c.height = crow;
  mat_c.width = ccol;
  mat_c.stride = ccol;

  //solve a*b
  clock_t begin = clock();
  MatMul(mat_a, mat_b, mat_c);
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;


  //output results
  FILE *f;
  f = fopen("distRes.txt", "w");
  for(i = 0; i<crow; i++){
    for(j = 0; j<ccol; j++){
      fprintf(f, "%lld ", mat_c.elements[i*ccol+j]);
    }
    fprintf(f, "\n");
  }
  fclose(f);

  printf("here\n");

  //free memory
  free(mat_a.elements);
  free(mat_b.elements);
  free(mat_c.elements);

  printf("\n--------------------------\nExecution took: %lf seconds\n", time_spent);

  return 0;
}
