#include <iostream>
using namespace std;

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
  cout << "What is the teperature of the left side of the rod for all time? ";
  cin >> tleft;
  cout << "What is the teperature of the right side of the rod for all time? ";
  cin >> tright;
  cout << "Initial conditions:" << endl;
  double oldt, newt[ndx];
  oldt[0] = tleft;
  oldt[ndx - 1] = tright;
  newt[0] = tleft;
  newt[ndx - 1] = tright;
  
}
