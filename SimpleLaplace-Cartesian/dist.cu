#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <algorithm>
using namespace std;
#define PI 3.14159265359

#define grid(i,j,ny) i*(ny+1)+j
#define omega 1.5
#define RelativeError 1e-4
#define epsilon 1e-10

#define threadx 16
#define thready 16

//-----------------------------------------------------------------------------
int getPos(int i, int j, int nx){
  return (((j-1) * (nx -1)) + i) - 1;
}

//-----------------------------------------------------------------------------
double getQ(int n, double temp, double lx, double ly){
  double res = cos(n*PI);
  res -= 1;
  res /= sinh(lx*n*PI/ly);
  res *= -2*temp;
  res /= (n*PI);
  return res;
}

//-----------------------------------------------------------------------------
double getT(double x, double y, double lx, double ly, double temp){
  double sum = 0;
  for(int n = 1; n<97; n++){
    sum += getQ(n, temp, lx, ly)*sinh(n*PI*x/ly)*sin(n*PI*y/ly);
  }
  return sum;
}

//-----------------------------------------------------------------------------
double findMaxTemp(double tleft, double tright, double ttop, double tbottom){
  double max = tleft;
  if(tright > max)
    max = tright;
  if(ttop > max)
    max = ttop;
  if(tbottom > max)
    max = tbottom;
  return max;
}

//-----------------------------------------------------------------------------
__global__ void update(double *grid_old, double * grid_new, bool isRed, bool *converge, int nx, int ny, double ldx, double ldy, int blocky){
  extern __shared__ double grid_s[];
  int i = blockIdx.x * threadx + threadIdx.x;//x position index
  int j = blockIdx.y * thready + threadIdx.y;//y position index
  int sharedPos = threadIdx.x * (thready+2) + threadIdx.y;//Making 2D into 1D for the shared memory position in this block
  int blockPos = 2*(blockIdx.x * blocky + blockIdx.y);
  if(isRed){//Could have just as well been !isRed
    blockPos++;//Because each block has two convergence flags, need to only update one of the two
  }
  /*if(blockPos%2 == 1 && (threadIdx.x == 0 || threadIdx.y == 0 || threadIdx.x == threadx+1 || threadIdx.y == thready+1)){
    printf("%d %d %d\n", blockPos, i, j);
  }*/
  if(i< nx + 1 && j < ny + 1){//Within the domain of the grid
    grid_s[sharedPos] = grid_old[grid(i,j,ny)];
  }
  __syncthreads();
  if(i< nx + 1 && j < ny + 1){//Within the domain of the grid
    converge[blockPos] = true;//Default. Then all you need is one 'false' to force another iteration
    if(i != 0 && i != nx && j != 0 && j!= ny){//boundaries
      if((i%2 == j%2) == isRed){//Red or not Red?
        if(threadIdx.x != 0 && threadIdx.y != 0 && threadIdx.x != threadx + 1 && threadIdx.y != thready + 1){//halo points
          grid_new[grid(i,j,ny)] = (1-omega)*grid_s[sharedPos];//copy a weighted fraction of the old
          //Then update with the remaining fraction with the new
          grid_new[grid(i,j,ny)] += omega *(
                                    grid_s[sharedPos-(thready+2)]*ldy*ldy/(2*(ldx*ldx+ldy*ldy)) + //left
                                    grid_s[sharedPos+thready+2]*ldy*ldy/(2*(ldx*ldx+ldy*ldy)) + //Right
                                    grid_s[sharedPos+1]*ldx*ldx/(2*(ldx*ldx+ldy*ldy)) + //top
                                    grid_s[sharedPos-1]*ldx*ldx/(2*(ldx*ldx+ldy*ldy)) //bottom
                                    );
          if((grid_new[grid(i,j,ny)] - grid_s[sharedPos])/(grid_s[sharedPos] + epsilon) > RelativeError){
            converge[blockPos] = false;
          }
        }//end of not halo
      }//end of red/black
    }//end of not boundaries
  }//end of domain of the grid
}//end of update

//-----------------------------------------------------------------------------
int main(){
  double tbot, ttop, tleft, tright;
  //cout << "Temperature of bottom side: ";
  cin >> tbot;
  //cout << "Temperature of top side: ";
  cin >> ttop;
  //cout << "Temperature of left side: ";
  cin >> tleft;
  //cout << "Temperature of right side: ";
  cin >> tright;
  int nx, ny;
  double lx, ly;
  //cout << "Number of segments in x: ";
  cin >> nx;
  //cout << "Number of segments in y: ";
  cin >> ny;
  //cout << "Length of x: ";
  cin >> lx;
  //cout << "Length of y: ";
  cin >> ly;
  double lnx = lx/nx;
  double lny = ly/ny;
  double *grid_h, *error, *percerror;
  double *grid_d_old, *grid_d_new;
  size_t gridSize = (nx + 1)*(ny + 1)*sizeof(double);
  cudaMalloc(&grid_d_new, gridSize);
  cudaMalloc(&grid_d_old, gridSize);
  grid_h = (double*)malloc(gridSize);
  error = (double*)malloc(gridSize);
  percerror = (double*)malloc(gridSize);

  //initialize grid with initial conditions
  //i is an index of x position
  for(int i = 0; i< nx + 1; i++){
    grid_h[i*(ny+1)] = tbot;
    grid_h[i*(ny+1)+ny] = ttop;
    //boundaries are defined so error = 0
    error[i*(ny+1)] = 0;
    error[i*(ny+1)+ny] = 0;
    percerror[i*(ny+1)] = 0;
    percerror[i*(ny+1)+ny] = 0;
  }

  for(int j = 0; j< ny + 1; j++){
    grid_h[j] = tleft;//i=0 here
    grid_h[nx*(ny+1)+j] = tright;//i=nx here
    //boundaries are defined so error = 0
    error[j] = 0;
    error[nx*(ny+1)+j] = 0;
    percerror[j] = 0;
    percerror[nx*(ny+1)+j] = 0;
  }

//Initial guess is based on a linear interpolation between the boundaries
  double deltai = (tright - tleft)/nx;
  double deltaj = (ttop - tbot)/ny;
  for(int i = 1; i< nx; i++){
    for(int j = 1; j < ny; j++){
        grid_h[i*(ny+1)+j] = 0.5*((tleft + deltai*i)+(tbot + deltaj*j));
    }
  }
  cudaMemcpy(grid_d_new, grid_h, gridSize, cudaMemcpyHostToDevice);

  dim3 threadDim(threadx + 2, thready + 2);
  int blockx = 1 + ((nx-1 -1)/threadx);//nx-1 is the number of interior x points
  int blocky = 1 + ((ny-1 - 1)/thready);//ny-1 is the number of interior y points
  dim3 blockDim(blockx, blocky);

  bool *converge_d;
  bool *converge_h;
  converge_h = (bool*)malloc(2*blockx*blocky*sizeof(bool));

  size_t convSize = 2*blockx*blocky*sizeof(bool);
  cudaMalloc(&converge_d, convSize);

  int steps = 0;
  bool didConverge = false;
  printf("Blocks: %d\n", blockx*blocky);
  while(!didConverge){
    //cout<<"BEFORE:"<<grid_h[7*(ny+1)+3]<<endl;
    //printf("Step: %d\n", steps);
    cudaMemcpy(grid_d_old, grid_d_new, gridSize, cudaMemcpyDeviceToDevice);
    //cudaDeviceSynchronize();
    update<<<blockDim,threadDim,(threadx+2)*(thready+2)*sizeof(double)>>>(grid_d_old, grid_d_new, true, converge_d, nx, ny, lnx, lny, blocky);
    //cudaDeviceSynchronize();
    cudaMemcpy(grid_d_old, grid_d_new, gridSize, cudaMemcpyDeviceToDevice);
    //cudaDeviceSynchronize();
    update<<<blockDim,threadDim,(threadx+2)*(thready+2)*sizeof(double)>>>(grid_d_old, grid_d_new, false, converge_d, nx, ny, lnx, lny, blocky);
    //cudaDeviceSynchronize();
    cudaMemcpy(converge_h, converge_d, convSize, cudaMemcpyDeviceToHost);
    //cudaDeviceSynchronize();
    steps++;

    didConverge = converge_h[0];
    for(int i = 1; i< 2*blockx*blocky; i++){
      didConverge = didConverge && converge_h[i];
    }
  }//converged
  printf("\nconverged in %d steps.\n", steps);
  cudaMemcpy(grid_h, grid_d_new, gridSize, cudaMemcpyDeviceToHost);
  //cudaDeviceSynchronize();
  //cout<<"AFTER:"<<grid_h[7*(ny+1)+3]<<endl;

  //fill grid
  for(int i = 1; i< nx; i++){
    for(int j = 1; j< ny; j++){
      double actual = getT(i*lnx, j*lny, lx, ly, findMaxTemp(tleft, tright, ttop, tbot));
      *(error+(i*(ny+1))+j) = *(grid_h+(i*(ny+1))+j) - actual;
      if(*(error+(i*(ny+1))+j) < 0){
        *(error+(i*(ny+1))+j) *= -1;
      }
      //cout << *(error+(i*(ny+1))+j) << " " << actual << " " << *(error+(i*(ny+1))+j)/actual << endl;
      *(percerror+(i*(ny+1))+j) = 100**(error+(i*(ny+1))+j)/(actual+epsilon);
    }
  }

/*
  cout << "\n-----------------" << endl;
  cout << "Final Results:" << endl;
  //see grid
  for(int j = ny; j> -1 ; j--){
    for(int i = 0; i< nx + 1; i++){
      cout << setprecision(3) << *(grid_h+(i*(ny+1))+j) << "\t";
    }
    cout << endl;
  }

  cout << "\n-----------------" << endl;
  cout << "Error:" << endl;
  //see grid
  for(int j = ny - 1; j> 0 ; j--){
    for(int i = 1; i< nx; i++){
      cout << setprecision(3) << *(error+(i*(ny+1))+j) << "\t";
    }
    cout << endl;
  }

  cout << "\n-----------------" << endl;
  cout << "Percent Error:" << endl;
  //see grid
  for(int j = ny - 1; j> 0 ; j--){
    for(int i = 1; i< nx; i++){
      cout << setprecision(3) << *(percerror+(i*(ny+1))+j) << "\t";
    }
    cout << endl;
  }*/
  /*//output data for plotly
  ofstream myfile;
  myfile.open("SimpleLaplaceCartesianDataforplotly2.txt");
  myfile << "Simple Laplace in Cartesian\n";
  myfile << "X (m)\n";
  myfile << "Y (m)\n";
  myfile << lnx << endl << lx << endl << lny << endl << ly << endl << nx << endl << ny << endl;
  for(int i = 0; i< nx + 1; i++){
    for(int j = 0; j<ny + 1; j++){
      myfile << *(grid+(i*(ny+1))+j) << ",";
    }
    myfile << endl;
  }

  myfile.close();
  */

  /*//output data for matlab
  ofstream myfile;
  myfile.open("matLabData.txt");
  for(int i = 0; i< nx + 1; i++){
    for(int j = 0; j<ny + 1; j++){
      myfile << i*lnx << " " << j*lny << " " << *(grid+(i*(ny+1))+j) << "\n";
    }
    myfile << endl;
  }
  myfile.close();
  */
  //Output data for contour maps plotly
  ofstream myfile;
  myfile.open("plotlyContour.txt");
  int i, j;
  for(i = 0; i< nx; i++){
    myfile << i*lnx << " ";
  }
  myfile << i*lnx << "\n";
  for(j = 0; j<ny; j++){
    myfile << j*lny << " ";
  }
  myfile << j*lny << "\n";

  for(j = 0; j<ny + 1; j++){
    for(i = 0; i< nx; i++){
      myfile << grid_h[i*(ny+1)+j] << " ";
    }
    myfile << grid_h[i*(ny+1)+j] << "\n";
  }
  myfile.close();

  cudaFree(grid_d_new);
  cudaFree(grid_d_old);
  cudaFree(converge_d);
  free(grid_h);
  free(error);
  free(percerror);
  free(converge_h);
}
