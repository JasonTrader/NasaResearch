#ifndef _GLOBALS_H_
#define _GLOBALS_H_

//NOTE globals should be declared in all caps to denote the fact that they are
//globals

#define PI 3.1415926535897932384//PI
#define MU0 4*PI*1e-7//mu naught
//NOTE in evals per block there might be more threads than the specidfied number for use as halo points
#define R_EVALS_PER_BLOCK 16//How many r values we wish to get results from in each block
#define Z_EVALS_PER_BLOCK 16//How many z values we wish to get results from in each block
#define EPSILON 1e-12//Prevents dividing by 0 when added to the denominator
#define q 1.60217662e-19//elementary charge
#define E0 8.854187817e-12//epsilon naught
#define kBoltz 1.38064852e-23//Boltzman's constant
#define Tp 20000//Temperature of positive ion
#define Tn 20000//Temperature of negative ion
#define Mp 1.6726219e-27//Mass of 1 AMU positive ion
#define Mn 1.6726219e-27//Mass of 1 AMU negative ion

//U is a vector of conserved variables
//This enum is placed in the order of which each of these conserved variables
//are placed in the U vector
enum U{
  massP,//Mass of positive ions at a particular space in the thruster
  massN,//Mass of negative ions at a particular space in the thruster
  momentumPR,//Momentum in r direction of positive ions at a particular space in the thruster
  momentumNR,//Momentum in r direction of negative ions at a particular space in the thruster
  momentumPZ,//Momentum in z direction of positive ions at a particular space in the thruster
  momentumNZ//Momentum in z direction of negative ions at a particular space in the thruster
};

//Old/new enum
enum on{
  o,//Old
  n//New
};

#endif
