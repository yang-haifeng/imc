#ifndef _GRID_H
#define _GRID_H

#include <math.h>
#include "typedef.h"
//#include "Eigen/Dense"

//using namespace Eigen;

class Grid{
  protected:
    int Nr, Ntheta;
    int NphiI, NthetaI;
    double * Density;
    double * BnuT;
    double * Stokes;

    double * rc, * rl, * rr;
    double * thetac, * thetal, * thetar;

    double * phiIc;
    double * thetaIc;
    double dphiI,dthetaI;

    void initDensity();
    void initBnuT();

    void iteration();
  public:
    Grid();
    ~Grid();
    double get_density(int ir, int it);
    double get_bnuT(int ir, int it);

    bool isInDomain(double r, double theta);
};

#endif
