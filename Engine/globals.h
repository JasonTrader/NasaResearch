#ifndef _GLOBALS_H_
#define _GLOBALS_H_

//NOTE globals should be declared in all caps to denote the fact that they are
//globals

#define PI 3.1415926535897932384
#define MU0 4*PI*1e-7
#define R_EVALS_PER_BLOCK 16
#define Z_EVALS_PER_BLOCK 16
#define EPSILON 1e-12

enum U{
  massP,
  massN,
  momentumPR,
  momentumNR,
  momentumPZ,
  momentumNZ
};

#endif
