#include <iostream>
using namespace std;

#define PI 3.14159
#define mu0 4*PI*0.0000001

int main(){
  double imax, rlength, eta, tstep, ldr;
  int numseg;
  cout << "What is your I max? ";
  cin >> imax;
  cout << "What is the length of your rod? ";
  cin >> rlength;
  cout << "What is eta? ";
  cin >> eta;
  cout << "How many segments would you like? ";
  cin >> numseg;
  ldr = rlength/numseg;
  cout << "What time step would you like? (must be less than " << 0.5*ldr*ldr*mu0/eta << " ) ";
  cin >> tstep;

  return 0;
}
