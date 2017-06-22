#ifndef _DUST_H
#define _DUST_H

#include <math.h>
#include <complex>
#include "typedef.h"

class Dust{ // Assuming small oblate dust grain for now.
  protected:
    std::complex<double> eps; // Complex dialectic constant
    double re; // Effective grain size in micron
    double s; // Axis ratio 
    std::complex<double> a1, a3;

    void calcAM();
  public:
    Dust();
    //Dust(std::complex<double> eps, double re, double s);
    ~Dust();
};

#endif
