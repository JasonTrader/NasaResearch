#include <iostream>
#include <iomanip>
using namespace std;


bool steadyState(double oldt[], double newt[], double thresh, int ndx){
  double totThresh = (ndx-2)*thresh;
  double sum = 0;
  for(int i = 1; i < ndx - 1; i++){
    double diff = newt[i] - oldt[i];
    if (diff < 0)
      diff *= -1;
    //sum += temp;
    if (diff > thresh)
      return(0);
  }
  return (1);
}

int main(){
  //lx is length of rod ndx is the number of segments and ldx is the length of each segment.
  double lx, ldx;
  int ndx;
  cout << "How long is your rod? ";
  cin >> lx;
  cout << "How many segments of your rod would you like? ";
  cin >> ndx;
  ldx = lx/ndx;
  double tleft, tright;
  cout << "What is the teperature of the left boundary of the rod initially? ";
  cin >> tleft;
  cout << "What is the teperature of the right boundary of the rod for all time? ";
  cin >> tright;
  cout << "what is your diffusivity constant? ";
  double difus;
  cin >> difus;
  cout << "What is you time step? (Must be smaller than " << (ldx*ldx)/(2*difus) << "s to ensure stability) " ;
  double tstep;
  cin >> tstep;
  cout << "What must the maximum change be to consider the solution as steady? ";
  double thresh;
  cin >> thresh;
  cout << "Initial conditions\n-----------------------" << endl;
  double oldt [ndx + 2];
  double newt[ndx + 2];
  oldt[0] = tleft;
  oldt[ndx + 1] = tright;
  newt[0] = tleft;
  newt[ndx + 1] = tright;

  for(int inc = 1; inc < ndx + 1; inc++){
    cout << "Segemt " << inc << ": ";
    cin >> newt[inc];
  }
  int tcount = 0;
  do{
    tcount++;
    oldt[0] = newt[0];
    for(int i = 1; i < (ndx + 1); i++){
      oldt[i] = newt[i];
      newt[i] = -2*oldt[i];
      newt[i] += oldt[i + 1];
      newt[i] += oldt[i - 1];
      newt[i] *= (difus*tstep/(ldx*ldx));
      newt[i] += oldt[i];
      if(i == 1)
        cout << newt[i] << " ";
      cout << newt[i] << " ";
    }
    newt[0] = newt[1];
    cout << newt[ndx + 1] << endl;
  } while(!steadyState(oldt, newt, thresh, ndx));
  cout << "\n---------------------------------\n";
  cout << "Solution converged in " << tcount * tstep << " seconds\n";
}
