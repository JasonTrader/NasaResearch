#ifndef _VOLTAGE_H_
#define _VOLTAGE_H_

#include <math.h>

double getInletVolt(double vmax, double vmin, double cycleTime, double currentTime){
  if(fmod(currentTime,cycleTime) < (cycleTime/2))
    return vmax;
  return vmin;
}


#endif
