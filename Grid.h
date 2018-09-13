#ifndef _GRID_H
#define _GRID_H

#include <iostream>
#include <math.h>
#include "typedef.h"
#include "utils.h"

#ifdef _MPI_
#include "mpi.h"
extern int world_size;
extern int my_rank;
#endif

class Grid{
  protected:
    int Nr, Ntheta;
    int NphiI, NthetaI;
    double * Density;
    double * BnuT;
    double * Stokes;

    double * Bfield;

    double * rc, * rl, * rr;
    double * thetac, * thetal, * thetar;

    double * phiIc;
    double * thetaIc;
    double dphiI,dthetaI;

    double r2min, r2max;
    double epsDS;

    double kappa_abs, kappa_sca, kappa_ext;

    void init();

    Vector calcScattering(int ir, int it, double x, double y, double z, 
    	double nx, double ny, double nz);
    Vector calcEmission(int ir, int it, double x, double y, double z, 
    	double nx, double ny, double nz);
    Matrix calcExtinction(int ir, int it, double x, double y, double z, 
    	double nx, double ny, double nz);

    Matrix muller_Matrix(double theta, double phi,
    	double nx, double ny, double nz, double phip, int ir, int it);

    void getSurface(double &x, double &y, double &z, double nx, double ny, 
        double nz, bool &status, int &ir, int &it);

    bool isInDomain(double x, double y, double z);

    void InitializeDensityWFunction(FieldFunction_ MyDensity);
    void InitializeBnuTWFunction(FieldFunction_ MyBnuT);

  public:
    Grid();
    ~Grid();
    double get_density(int ir, int it);
    double get_bnuT(int ir, int it);

    void iteration(bool ScaFlag=true);

    void saveStokes(std::string fName="stokes.bin");
    void loadStokes(std::string fName="stokes.bin");

    void Image(double inc, int Npix=100, bool ifsca=true, std::string fName="image.out");
    Vector OnePointImage(double x0, double y0, double inc);
    Vector OnePointImage_wID(double inc, int ID, int Npix=100, bool ifsca=true);
    Vector Integrate(double x0,double y0,double z0, 
        double nx,double ny,double nz, bool ScaFlag=true);
    void moveOneCell(double x, double y, double z, double nx, double ny, double nz,
    	double &ds, int &ir, int &it);
    bool findCell(double x, double y, double z, int &ir, int &it);
};

#endif
