#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define PI 3.1415926535897932384
#define mu0 4*PI*1e-7

int main(){
  double Il, Ir, rlength, eta, tstep, ldr, ldz, tottime, zlength;
  int rseg, zseg;
  printf("What is your left I? ");
  scanf("%lf", &Il);
  printf("What is your right I? ");
  scanf("%lf", &Ir);
  printf("What is the radius of your rod? ");
  scanf("%lf", &rlength);
  printf("What is the length of your rod? ");
  scanf("%lf", &zlength);
  printf("What is eta? ");
  scanf("%lf", &eta);
  printf("How many segments would you like per radius? ");
  scanf("%d", &rseg);
  printf("How many segments would you like per length? ");
  scanf("%d", &zseg);
  ldr = rlength/(rseg+1);
  ldz = zlength/(zseg+1);
  double smallest = ldr;
  if(ldz < ldr)
    smallest = ldz;
  tstep = 0.125*smallest*smallest*mu0/eta;
  printf("How long would you like to run? ");
  scanf("%lf", &tottime);

  double dI = (Ir - Il) / (zseg+2);

  double *grid_new, *grid_old;
  size_t grid_size = (rseg+2)*(zseg+2)*sizeof(double);
  grid_new = (double*)malloc(grid_size);
  grid_old = (double*)malloc(grid_size);

  long int i;

  for(i = 0; i<(rseg+2)*(zseg+2); i++){
    long int rrow = i/(zseg+2);
    long int zcol = i%(zseg+2);
    grid_new[i] = (1-(rrow*rrow*ldr*ldr/(3*rlength*rlength)))*3*mu0*(Il + dI*zcol)*rrow*ldr/(4*PI*rlength*rlength);
  }

  FILE *myfile;
  myfile = fopen("init.txt", "w");



  for(i = 0; i< zseg+1; i++)
    fprintf(myfile, "%lf ", i*ldz);
  fprintf(myfile, "%lf\n", i*ldz);

  for(i = 0; i< rseg+1; i++)
    fprintf(myfile, "%lf ", i*ldr);
  fprintf(myfile, "%lf\n", i*ldr);

  for(i = 0; i< (rseg + 2)*(zseg+2); i++){
    if(i%(zseg+2)==zseg+1)
      fprintf(myfile, "%lf\n", grid_new[i]);
    else
      fprintf(myfile, "%lf ", grid_new[i]);

  }

  fclose(myfile);

  double r_aug = eta*tstep/(mu0*ldr*ldr);
  double z_aug = eta*tstep/(mu0*ldz*ldz);
  long int total_steps = tottime / tstep;
  printf("\nSteps: %ld\n", total_steps);


  clock_t begin, end;
  double time_spent;
  begin = clock();
  long int tcount = 0;

  while((tcount*tstep) < tottime){
    tcount++;

    //copy new to old
    for(i = 0; i<(rseg + 2)*(zseg+2); i++){
      grid_old[i] = grid_new[i];
    }

    //update
    for(i = 0; i<(rseg + 2)*(zseg+2); i++){
      int r = i/(zseg+2);
      int z = i%(zseg+2);

      if(z != 0 && z != zseg+1){
        if(r!= 0 && r!= rseg+1){
          if(r==1){
            grid_new[i] += r_aug*(2*grid_old[i+(zseg+2)] - 4*grid_old[i]) +
              z_aug * (grid_old[i+1] - 2*grid_old[i] + grid_old[i-1]);
          }
          else{
            grid_new[i] += r_aug*((1+(1/(2*i)))*grid_old[i+(zseg+2)] + (-2-(1/(i*i)))*grid_old[i] + (1-(1/(2*i)))*grid_old[i-(zseg+2)])
              +z_aug*(grid_old[i+1] - 2*grid_old[i] + grid_old[i-1]);
          }
        }
      }
    }
  }
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;


  myfile = fopen("res.txt", "w");

  for(i = 0; i< zseg+1; i++)
    fprintf(myfile, "%lf ", i*ldz);
  fprintf(myfile, "%lf\n", i*ldz);

  for(i = 0; i< rseg+1; i++)
    fprintf(myfile, "%lf ", i*ldr);
  fprintf(myfile, "%lf\n", i*ldr);

  for(i = 0; i< (rseg + 2)*(zseg+2); i++){
    if(i%(zseg+2)==zseg+1)
      fprintf(myfile, "%lf\n", grid_new[i]);
    else
      fprintf(myfile, "%lf ", grid_new[i]);

  }

  fclose(myfile);

  free(grid_new);
  free(grid_old);



  printf("\n------------------------------------\nExecution took: %lf sec\n", time_spent);


  return 0;
}
