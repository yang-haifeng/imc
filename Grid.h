#ifndef _GRID_H
#define _GRID_H

#include <iostream>
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
    double * Stokes1;

    double * rc, * rl, * rr;
    double * thetac, * thetal, * thetar;

    double * phiIc;
    double * thetaIc;
    double dphiI,dthetaI;

    double r2max;

    double kappa_abs, kappa_sca, kappa_ext;

    void initDensity();
    void initBnuT();

    void moveOneCell(double x, double y, double z, double nx, double ny, double nz,
    	double &ds, int &ir, int &it);
    void calc_Scattering(int ir, int it, double x, double y, double z, 
    	double nx, double ny, double nz, double tau, double dtau,
	double &I, double &Q, double &U, double &V);

    void muller_Matrix(double theta, double phi, double nx, double ny, double nz,
    	double I, double Q, double U, double V,
    	double &dI, double &dQ, double &dU, double &dV,
	int ir, int it);

  public:
    Grid();
    ~Grid();
    double get_density(int ir, int it);
    double get_bnuT(int ir, int it);

    bool isInDomain(double x, double y, double z);

    void zeroIter();
    void iteration();

    void saveStokes(std::string fName="stokes.bin");
};

#endif
