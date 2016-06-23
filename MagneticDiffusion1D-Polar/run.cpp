#include <iostream>
#include <stdio.h>
#include <time.h>
#include <iomanip>
using namespace std;

#define PI 3.1415926535897932384
#define mu0 4*PI*1e-7

bool converge(double rod_new[], double rod_old[], double length, double thresh){
  for(int i=1; i<length+1; i++){
    double diff = rod_old[i] - rod_new[i];
    if(diff < 0)
      diff*= -1;
    if(diff > thresh)
      return false;
  }
  return true;
}

int main(){

  clock_t begin, end;
  double time_spent;
  begin = clock();

  FILE *myfile;
  myfile = fopen("results.txt", "w");
  double imax, rlength, eta, tstep, ldr, conv;
  int numseg;
  cout << "What is your I max? ";
  cin >> imax;
  cout << "What is the length of your rod? ";
  cin >> rlength;
  cout << "What is eta? ";
  cin >> eta;
  cout << "How many segments would you like? ";
  cin >> numseg;
  ldr = rlength/(numseg+1);
  cout << "What time step would you like? (must be less than " << 0.5*ldr*ldr*mu0/eta << " ) ";
  cin >> tstep;
  cout << "What is the threshold for your convergence? ";
  cin >> conv;

  double rod_new[numseg+2];
  double rod_old[numseg+2];
  rod_new[0] = 0;
  rod_new[numseg + 1] = mu0*imax/(2*PI*rlength);
  rod_old[0] = 0;
  rod_old[numseg + 1] = mu0*imax/(2*PI*rlength);
  for(int i = 1; i<numseg+1; i++){
    rod_new[i] = (1-(i*i*ldr*ldr/(3*rlength*rlength)))*3*mu0*imax*i*ldr/(4*PI*rlength*rlength);
  }

  int out;
  //output r values
  for(out = 0; out<numseg+1; out++){
    fprintf( myfile, "%lf ", out*ldr );
}
fprintf( myfile, "%lf\n", out*ldr );

  double aug = eta*tstep/(mu0*ldr*ldr);
  int tcount = 0;
  do{
    if(tcount%10000==0){
      for(out = 0; out<numseg+1; out++){
        fprintf( myfile, "%lf ", *(rod_new+out) );
    }
    fprintf( myfile, "%lf\n", *(rod_new+out) );

    }
    tcount++;

    //copy new to old
    for(int i = 1; i<numseg+1; i++){
      rod_old[i] = rod_new[i];
    }

    //update

    rod_new[1]+= aug*(2*rod_old[2] - 4*rod_old[1]);
    for(int i = 2; i<numseg+1; i++){
      rod_new[i] += aug*((1+(1/(2*i)))*rod_old[i+1] + (-2-(1/(i*i)))*rod_old[i] + (1-(1/(2*i)))*rod_old[i-1]);
    }
  } while(!converge(rod_new, rod_old, numseg, conv));

  fprintf(myfile, "STOP\n");
  fclose(myfile);

  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  cout << "\n------------------------------------\nExecution took: "<<  time_spent << " sec\n";

  return 0;
}
