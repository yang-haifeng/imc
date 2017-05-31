#include "Grid.h"

#define Ntot 50

Grid::Grid(){
  Nr = Ntot; Ntheta = Ntot; // Number in spacial grid.
  NphiI = Ntot; NthetaI = Ntot; // Number in scattering angle grid.

  // Allocate memory for the arrays.
  // Density and BnuT are 2D array. ir*Ntheta+itheta
  Density = new double [Nr*Ntheta];
  BnuT = new double [Nr*Ntheta];
  // Stokes number is a 5D array.
  // ir*Ntheta*NphiI*NthetaI*4 + itheta*NphiI*NthetaI*4 
  // + iphiI*NthetaI*4 + ithetaI*4 + i(0,1,2,3)
  // i(0,1,2,3) for Stokes I, Q, U, V
  Stokes = new double [Nr*Ntheta*NphiI*NthetaI*4];
  // Radial grid. 1D arrays for center, left and right edge of the cell;
  rc = new double [Nr];
  rl = new double [Nr];
  rr = new double [Nr];

  // Longitudual grid. 1D arrays for center, left and right edge of the cell;
  thetac = new double [Ntheta];
  thetal = new double [Ntheta];
  thetar = new double [Ntheta];

  // Angular grid.
  phiIc = new double [NphiI];
  thetaIc = new double [NthetaI];

  this->initDensity();
  this->initBnuT();

  for(int i=0;i<Nr*Ntheta*NphiI*NthetaI*4; i++) Stokes[i]=0;

  for(int i=0;i<Nr;i++){
    rc[i] = i*10*AU+5*AU;
    rl[i] = i*10*AU;
    rr[i] = i*10*AU+10.*AU;
  }
  for(int i=0;i<Ntheta;i++){ // Max theta is PI now. Note that it is pretty bad.
    // I'll work on a modification later. 
    thetac[i] = PI/Ntheta*(i+0.5);
    thetal[i] = PI/Ntheta*(i);
    thetar[i] = PI/Ntheta*(i+1.);
  }

  // Angular grid and its step sizes.
  for(int i=0;i<NphiI;i++) phiIc[i] = 2.*PI/NphiI*(i+0.5);
  dphiI = 2.*PI/NphiI;
  for(int i=0;i<NthetaI;i++) thetaIc[i] = PI/2./NthetaI*(i+0.5);
  dthetaI = PI/2./NthetaI;

  r2max = rr[Nr-1]*rr[Nr-1];

  kappa_abs = 1.;
  kappa_sca = 0.1;
  kappa_ext = kappa_abs+kappa_sca;
}

Grid::~Grid(){
  delete [] Density;
  delete [] BnuT;
  delete [] Stokes;

  delete [] rc;
  delete [] rl;
  delete [] rr;

  delete [] thetac;
  delete [] thetal;
  delete [] thetar;

  delete [] phiIc;
  delete [] thetaIc;
}
