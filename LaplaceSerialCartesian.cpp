#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;
#define PI 3.14159265359

int getPos(int i, int j, int nx){
  return (((j-1) * (nx -1)) + i) - 1;
}

double getQ(int n, double temp, double lx, double ly){
  /*double res = cos(PI*n);
  res *= ratio;
  res = ratio - res;
  res /= sinh(PI*n);
  res *= 2;
  res /= n;
  */
  double res = cos(n*PI);
  res -= 1;
  res /= sinh(lx*n*PI/ly);
  res *= -2*temp;
  res /= (n*PI);
  return res;
}

double getT(double x, double y, double lx, double ly, double temp){
  double sum = 0;
  for(int n = 1; n<97; n++){
    sum += getQ(n, temp, lx, ly)*sinh(n*PI*x/ly)*sin(n*PI*y/ly);
  }
  return sum;
}

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

int main(){
  double tbot, ttop, tleft, tright;
  cout << "Temperature of bottom side: ";
  cin >> tbot;
  cout << "Temperature of top side: ";
  cin >> ttop;
  cout << "Temperature of left side: ";
  cin >> tleft;
  cout << "Temperature of right side: ";
  cin >> tright;
  int nx, ny;
  double lx, ly;
  cout << "Number of segments in x: ";
  cin >> nx;
  cout << "Number of segments in y: ";
  cin >> ny;
  cout << "Length of x: ";
  cin >> lx;
  cout << "Length of y: ";
  cin >> ly;
  double lnx = lx/nx;
  double lny = ly/ny;
  double grid[nx + 1][ny + 1];
  double error[nx + 1][ny + 1];
  double percerror[nx + 1][ny + 1];
  //initialize grid with initial conditions
  for(int i = 0; i< nx + 1; i++){
    grid[i][0] = tbot;
    grid[i][ny] = ttop;
    error[i][0] = 0;
    error[i][ny] = 0;
    percerror[i][0] = 0;
    percerror[i][ny] = 0;
  }
  for(int j = 0; j< ny + 1; j++){
    grid[nx][j] = tright;
    grid[0][j] = tleft;
    error[nx][j] = 0;
    error[0][j] = 0;
    percerror[nx][j] = 0;
    percerror[0][j] = 0;
  }


  int matNum = nx - 1;
  matNum *= (ny -1);
  double mat[matNum][matNum + 1];

  int counter = matNum;
  double aug = lnx / lny;
  aug *= aug;

  while(counter --> 0){
    int j = (counter / (nx-1)) + 1;
    int i = (counter % (nx-1)) + 1;
    //zero out matrix
    for(int q = 0; q< (matNum + 1); q++){
      mat[counter][q] = 0;
    }

    //fill matrix
    //matPos = ((j-1) * (nx -1)) + i;
    int matPos = counter;
    if(i == (nx - 1)){
      mat[counter][matNum] -= grid[nx][j];
    }
    else{
      mat[counter][matPos + 1] = 1;
    }
    if(i == 1){
      mat[counter][matNum] -= grid[0][j];
    }
    else{
      mat[counter][matPos - 1] = 1;
    }
    if(j == (ny - 1)){
      mat[counter][matNum] -= (aug * grid[i][ny]);
    }
    else{
      mat[counter][matPos + (nx - 1)] = aug;
    }
    if(j == 1){
      mat[counter][matNum] -= (aug * grid[i][0]);
    }
    else{
      mat[counter][matPos - (nx - 1)] = aug;
    }
    mat[counter][counter] = -1*((2*aug) + 2);
  }
  //solve matrix
  //forward
  for(int q = 0; q<matNum; q++){
    if(mat[q][q] != 1){
      double temp = mat[q][q];
      for(int w = 0; w < matNum + 1; w++){
        mat[q][w] /= temp;
      }
    }
    for(int incq = q + 1; incq < matNum; incq++){
      double factor = -1*(mat[incq][q]/mat[q][q]);
      for(int incw = q; incw < matNum + 1; incw++){
        mat[incq][incw] += (factor * mat[q][incw]);
      }
    }
  }

  //backward
  for(int q = matNum - 1; q > -1; q--){
    for(int decq = q - 1; decq > -1; decq--){
      double factor = -1*(mat[decq][q]/mat[q][q]);
      mat[decq][q] = 0;
      mat[decq][matNum] += factor * mat[q][matNum];
    }
  }

  //fill grid
  for(int i = 1; i< nx; i++){
    for(int j = 1; j< ny; j++){
      grid[i][j] = mat[getPos(i, j, nx)][matNum];
      double actual = getT(i*lnx, j*lny, lx, ly, findMaxTemp(tleft, tright, ttop, tbot));
      error[i][j] = grid[i][j] - actual;
      if(error[i][j] < 0){
        error[i][j] *= -1;
      }
      percerror[i][j] = 100*error[i][j]/actual;
    }
  }

  /*//see matrix
  for(int q = 0; q< matNum; q++){
    for(int w = 0; w< matNum + 1; w++){
      cout << mat[q][w] << " ";
    }
    cout << endl;
  }*/



  cout << "\n-----------------" << endl;
  cout << "Final Results:" << endl;
  //see grid
  for(int j = ny; j> -1 ; j--){
    for(int i = 0; i< nx + 1; i++){
      cout << setprecision(3) << grid[i][j] << "\t";
    }
    cout << endl;
  }

  cout << "\n-----------------" << endl;
  cout << "Error:" << endl;
  //see grid
  for(int j = ny; j> -1 ; j--){
    for(int i = 0; i< nx + 1; i++){
      cout << setprecision(3) << error[i][j] << "\t";
    }
    cout << endl;
  }

  cout << "\n-----------------" << endl;
  cout << "Percent Error:" << endl;
  //see grid
  for(int j = ny; j> -1 ; j--){
    for(int i = 0; i< nx + 1; i++){
      cout << setprecision(2) << percerror[i][j] << "\t";
    }
    cout << endl;
  }



}
