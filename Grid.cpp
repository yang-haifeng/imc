#include "Grid.h"

#define Ntot 10
#define RMAX 100

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
  // Bfield is 3D array. The same 2D as Density + 1D of length 3.
  Bfield = new double [Nr*Ntheta*3];
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

  for(int i=0;i<Nr*Ntheta*NphiI*NthetaI*4; i++)
    Stokes[i]=0; 

  double dr = RMAX*AU/Nr;
  for(int i=0;i<Nr;i++){
    rc[i] = dr*(i+0.5);
    rl[i] = dr*i;
    rr[i] = dr*(i+1);
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
  for(int i=0;i<NthetaI;i++) thetaIc[i] = PI/NthetaI*(i+0.5);
  dthetaI = PI/NthetaI;

  r2min = rl[0]*rl[0];
  r2max = rr[Nr-1]*rr[Nr-1];

  kappa_abs = 1.;
  kappa_sca = 0.1;
  kappa_ext = kappa_abs+kappa_sca;

  this->initDensity();
  this->initBnuT();
  this->initBfield();
}

Grid::~Grid(){
  delete [] Density;
  delete [] BnuT;
  delete [] Stokes;

  delete [] Bfield;

  delete [] rc;
  delete [] rl;
  delete [] rr;

  delete [] thetac;
  delete [] thetal;
  delete [] thetar;

  delete [] phiIc;
  delete [] thetaIc;
}
