#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;
#define PI 3.14159265359

int getPos(int i, int j, int nr){
  return (((j-1) * (nr -1)) + i) - 1;
}

double getQ(int n, double temp, double lr, double lz){
  double res = cos(n*PI);
  res -= 1;
  res /= sinh(lr*n*PI/lz);
  res *= -2*temp;
  res /= (n*PI);
  return res;
}

double getT(double x, double y, double lr, double lz, double temp){
  double sum = 0;
  for(int n = 1; n<97; n++){
    sum += getQ(n, temp, lr, lz)*sinh(n*PI*x/lz)*sin(n*PI*y/lz);
  }
  return sum;
}

double findMaxTemp(double tleft, double tright, double tout){
  double max = tleft;
  if(tright > max)
    max = tright;
  if(tout > max)
    max = tout;
  return max;
}

int main(){
  double tout, tleft, tright;
  cout << "Temperature of outside: ";
  cin >> tout;
  cout << "Temperature of left side: ";
  cin >> tleft;
  cout << "Temperature of right side: ";
  cin >> tright;
  int nr, nz;
  double lr, lz;
  cout << "Number of segments in x: ";
  cin >> nr;
  cout << "Number of segments in y: ";
  cin >> nz;
  cout << "Length of x: ";
  cin >> lr;
  cout << "Length of y: ";
  cin >> lz;
  double lnr = lr/nr;
  double lnz = lz/nz;
  double grid[nr + 1][nz + 1];
  double error[nr + 1][nz + 1];
  double percerror[nr + 1][nz + 1];
  //initialize grid with initial conditions
  for(int i = 0; i< nr + 1; i++){
    grid[i][0] = tleft;
    grid[i][nz] = tright;
    error[i][0] = 0;
    error[i][nz] = 0;
    percerror[i][0] = 0;
    percerror[i][nz] = 0;
  }
  for(int j = 0; j< nz + 1; j++){
    grid[nr][j] = tout;
    error[nr][j] = 0;
    error[0][j] = 0;
    percerror[nr][j] = 0;
    percerror[0][j] = 0;
  }


  int matNum = nr - 1;
  matNum *= (nz -1);
  double mat[matNum][matNum + 1];

  int counter = matNum;
  double vert_aug = lnr*lnr/(lnz*lnz);

  while(counter --> 0){
    int j = (counter / (nr-1)) + 1;
    int i = (counter % (nr-1)) + 1;
    //zero out matrix
    for(int q = 0; q< (matNum + 1); q++){
      mat[counter][q] = 0;
    }

    //fill matrix
    if(i == (nr - 1)){
      mat[counter][matNum] -= (1+(1/(2*i)))*grid[nr][j];
    }
    else{
      mat[counter][counter + 1] = 1+(1/(2*i));
    }
    if(i == 1){
      mat[counter][counter] += 0.5;
    }
    else{
      mat[counter][counter - 1] = 1-(1/(2*i));
    }
    if(j == (nz - 1)){
      mat[counter][matNum] -= (vert_aug * grid[i][nz]);
    }
    else{
      mat[counter][counter + (nr - 1)] = vert_aug;
    }
    if(j == 1){
      mat[counter][matNum] -= (vert_aug * grid[i][0]);
    }
    else{
      mat[counter][counter - (nr - 1)] = vert_aug;
    }
    mat[counter][counter] += -1*((2*vert_aug) + 2);
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
  for(int i = 1; i< nr; i++){
    for(int j = 1; j< nz; j++){
      grid[i][j] = mat[getPos(i, j, nr)][matNum];
      double actual = getT(i*lnr, j*lnz, lr, lz, findMaxTemp(tleft, tright, tout));
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
  for(int i = 0; i< nr + 1; i++){
    for(int j = 0; j < nr + 1 ; j++){
      cout << setprecision(3) << grid[i][j] << "\t";
    }
    cout << endl;
  }

  cout << "\n-----------------" << endl;
  cout << "Error:" << endl;
  //see grid
  for(int i = 0; i< nr + 1; i++){
    for(int j = 0; j < nr + 1 ; j++){
      cout << setprecision(3) << error[i][j] << "\t";
    }
    cout << endl;
  }

  cout << "\n-----------------" << endl;
  cout << "Percent Error:" << endl;
  //see grid
  for(int i = 0; i< nr + 1; i++){
    for(int j = 0; j < nr + 1 ; j++){
      cout << setprecision(2) << percerror[i][j] << "\t";
    }
    cout << endl;
  }



}
