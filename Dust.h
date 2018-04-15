#ifndef _DUST_H
#define _DUST_H

#include <math.h>
#include <complex>
#include "typedef.h"
//#include "Grid.h"

class Grid;

class Dust{ // Assuming small oblate dust grain for now.
    friend class Grid;
  protected:
    std::complex<double> eps; // Complex dialectic constant
    double re; // Effective grain size in micron
    double s; // Axis ratio 
    std::complex<double> a1, a3;
    double p0;
  public:
    Dust();
    Dust(std::complex<double>, double, double);
    ~Dust();
    void calcAM();
};

#endif
