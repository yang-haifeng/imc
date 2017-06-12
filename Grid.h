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

    double * Bfield;

    double * rc, * rl, * rr;
    double * thetac, * thetal, * thetar;

    double * phiIc;
    double * thetaIc;
    double dphiI,dthetaI;

    double r2min, r2max;

    double kappa_abs, kappa_sca, kappa_ext;

    void initDensity();
    void initBnuT();
    void initBfield();

    void Integrate(double x0, double y0, double z0, double nx, double ny, double nz,
        double &I, double &Q, double &U, double &V, bool ScaFlag=true);
    void moveOneCell(double x, double y, double z, double nx, double ny, double nz,
    	double &ds, int &ir, int &it);
    void calc_Scattering(int ir, int it, double x, double y, double z, 
    	double nx, double ny, double nz, double tau, double dtau,
	double &I, double &Q, double &U, double &V);
    void calcEmission(int ir, int it, double x, double y, double z, 
    	double nx, double ny, double nz, 
	double &dI, double &dQ, double &dU, double &dV);

    void muller_Matrix(double theta, double phi, double nx, double ny, double nz,
    	double I, double Q, double U, double V, double phip,
    	double &dI, double &dQ, double &dU, double &dV,
	int ir, int it);

    void getSurface(double &x, double &y, double &z, double nx, double ny, 
        double nz, bool &status, int &ir, int &it);

    bool isInDomain(double x, double y, double z);
    bool findCell(double x, double y, double z, int &ir, int &it);

  public:
    Grid();
    ~Grid();
    double get_density(int ir, int it);
    double get_bnuT(int ir, int it);

    void iteration(bool ScaFlag=true);

    void saveStokes(std::string fName="stokes.bin");

    void Image(double inc, int Npix=100, std::string fName="image.out");
};

#endif
