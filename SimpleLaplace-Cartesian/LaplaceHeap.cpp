#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;
#define PI 3.14159265359

int getPos(int i, int j, int nx){
  return (((j-1) * (nx -1)) + i) - 1;
}

double getQ(int n, double temp, double lx, double ly){
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
  double *grid = new double[(nx + 1)*(ny + 1)];
  double *error = new double[(nx + 1)*(ny + 1)];
  double *percerror = new double[(nx + 1)*(ny + 1)];
  //initialize grid with initial conditions
  for(int i = 0; i< (nx + 1)*(ny+1); i++){
    if((i%(ny+1)) == 0){
      *(grid+i) = tbot;
      *(error+i) = 0;
      *(percerror+i) = 0;
    }
    if((i/(ny+1)) == 0){
      *(grid+i) = tleft;
      *(error+i) = 0;
      *(percerror+i) = 0;
    }
    if(i%(ny+1) == ny){
      *(grid+i) = ttop;
      *(error+i) = 0;
      *(percerror+i) = 0;
    }
    if(i/(ny+1) == nx){
      *(grid+i) = tright;
      *(error+i) = 0;
      *(percerror+i) = 0;
    }
  }

  int matNum = nx - 1;
  matNum *= (ny -1);
  double *mat = new double[matNum*(matNum + 1)];

  int counter = matNum;
  double aug = lnx / lny;
  aug *= aug;

  while(counter --> 0){
    int j = (counter / (nx-1)) + 1;
    int i = (counter % (nx-1)) + 1;
    //zero out matrix
    for(int q = 0; q< (matNum + 1); q++){
      *(mat+(counter*(matNum + 1))+q) = 0;
    }

    //fill matrix
    //counter = ((j-1) * (nx -1)) + i;
    if(i == (nx - 1)){
      *(mat+(counter*(matNum + 1))+matNum) -= *(grid+(nx*(ny+1))+j);
    }
    else{
      *(mat+(counter*(matNum + 1))+(counter + 1)) = 1;
    }
    if(i == 1){
      *(mat+(counter*(matNum + 1))+matNum) -= *(grid+j);
    }
    else{
      *(mat+(counter*(matNum + 1))+(counter - 1)) = 1;
    }
    if(j == (ny - 1)){
      *(mat+(counter*(matNum + 1))+matNum) -= (aug * *(grid+(i*(ny+1))+ny));
    }
    else{
      *(mat+(counter*(matNum + 1))+(counter + (nx - 1))) = aug;
    }
    if(j == 1){
      *(mat+(counter*(matNum + 1))+matNum) -= aug * (*(grid+(i*(ny+1))));
    }
    else{
      *(mat+(counter*(matNum + 1))+(counter - (nx - 1))) = aug;
    }
    *(mat+(counter*(matNum + 1))+counter) = -1*((2*aug) + 2);
  }

  //solve matrix
  //forward
  for(int q = 0; q<matNum; q++){
    if(*(mat+(q*(matNum + 1))+q) != 1){
      double temp = *(mat+(q*(matNum + 1))+q);
      for(int w = q; w < matNum + 1; w++){
        *(mat+(q*(matNum + 1))+w) /= temp;
      }
    }
    for(int incq = q + 1; incq < matNum; incq++){
      double factor = -1*(*(mat+(incq*(matNum+1))+q));
      for(int incw = q; incw < matNum + 1; incw++){
        *(mat+(incq*(matNum+1))+incw) += factor * (*(mat+(q*(matNum+1))+incw));
      }
    }
  }

  //backward
  for(int q = matNum - 1; q > -1; q--){
    for(int decq = q - 1; decq > -1; decq--){
      double factor = -1*(*(mat+(decq*(matNum+1))+q));
      *(mat+(decq*(matNum+1))+q) = 0;
      *(mat+(decq*(matNum+1))+matNum) += factor * (*(mat+(q*(matNum+1))+matNum));
    }
  }

  //fill grid
  for(int i = 1; i< nx; i++){
    for(int j = 1; j< ny; j++){
      *(grid+(i*(ny+1))+j) = *(mat+(getPos(i, j, nx)*(matNum+1))+matNum);
      double actual = getT(i*lnx, j*lny, lx, ly, findMaxTemp(tleft, tright, ttop, tbot));
      *(error+(i*(ny+1))+j) = *(grid+(i*(ny+1))+j) - actual;
      if(*(error+(i*(ny+1))+j) < 0){
        *(error+(i*(ny+1))+j) *= -1;
      }
      *(percerror+(i*(ny+1))+j) = 100*(*(grid+(i*(ny+1))+j))/actual;
    }
  }

  //see matrix
  //for(int q = 0; q< matNum; q++){
    //for(int w = 0; w< matNum + 1; w++){
      //cout << *(mat+(q*(matNum+1))+w) << " ";
    //}
    //cout << endl;
  //}


/*
  cout << "\n-----------------" << endl;
  cout << "Final Results:" << endl;
  //see grid
  for(int j = ny; j> -1 ; j--){
    for(int i = 0; i< nx + 1; i++){
      cout << setprecision(3) << *(grid+(i*(ny+1))+j) << "\t";
    }
    cout << endl;
  }

  cout << "\n-----------------" << endl;
  cout << "Error:" << endl;
  //see grid
  for(int j = ny; j> -1 ; j--){
    for(int i = 0; i< nx + 1; i++){
      cout << setprecision(3) << *(error+(i*(ny+1))+j) << "\t";
    }
    cout << endl;
  }

  cout << "\n-----------------" << endl;
  cout << "Percent Error:" << endl;
  //see grid
  for(int j = ny; j> -1 ; j--){
    for(int i = 0; i< nx + 1; i++){
      cout << setprecision(2) << *(grid+(i*(ny+1))+j) << "\t";
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

  //output data for matlab
  ofstream myfile;
  myfile.open("matLabData.txt");
  for(int i = 0; i< nx + 1; i++){
    for(int j = 0; j<ny + 1; j++){
      myfile << i*lnx << " " << j*lny << " " << *(grid+(i*(ny+1))+j) << "\n";
    }
    myfile << endl;
  }

  myfile.close();
}
