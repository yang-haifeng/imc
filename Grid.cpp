#include "Grid.h"

Grid::Grid(){
  Nr = 10; Ntheta = 10; // Number in spacial grid.
  NphiI = 10; NthetaI = 10; // Number in scattering angle grid.

  // Allocate memory for the arrays.
  // Density and BnuT are 2D array. ir*Ntheta+itheta
  Density = new double [Nr*Ntheta];
  BnuT = new double [Nr*Ntheta];
  // Stokes number is a 5D array.
  // ir*Ntheta*NphiI*NthetaI*4 + itheta*NphiI*NthetaI*4 
  // + iphiI*NthetaI*4 + ithetaI*4 + i(0,1,2,3)
  // i(0,1,2,3) for Stokes I, Q, U, V
  Stokes = new double [Nr*Ntheta*NphiI*NthetaI*4];
  // Center of the grid. All four arrays are 1D.
  rc = new double [Nr];
  thetac = new double [Ntheta];
  phiIc = new double [NphiI];
  thetaIc = new double [NthetaI];
  // Left and right side of the grid.
  rl = new double [Nr];
  rr = new double [Nr];
  thetal = new double [Ntheta];
  thetar = new double [Ntheta];

  this->initDensity();
  this->initBnuT();

  for(int i=0;i<Nr*Ntheta*NphiI*NthetaI*4; i++) Stokes[i]=0;

  for(int i=0;i<Nr;i++){
    rc[i] = i*10*AU+5*AU;
    rl[i] = i*10*AU;
    rr[i] = i*10*AU+10.*AU;
  }
  for(int i=0;i<Ntheta;i++){
    thetac[i] = PI/2./Ntheta*(i+0.5);
    thetal[i] = PI/2./Ntheta*(i);
    thetar[i] = PI/2./Ntheta*(i+1.);
  }

  for(int i=0;i<NphiI;i++) phiIc[i] = 2.*PI/Ntheta*(i+0.5);
  for(int i=0;i<NthetaI;i++) thetaIc[i] = PI/2./Ntheta*(i+0.5);

  // Step size in angular grid.

}

Grid::~Grid(){
  delete [] Density;
  delete [] BnuT;
  delete [] Stokes;
  delete [] rc;
  delete [] thetac;
  delete [] rl;
  delete [] thetal;
  delete [] rr;
  delete [] thetar;
  delete [] phiIc;
  delete [] thetaIc;
}
