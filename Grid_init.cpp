#include "Grid.h"

#define RMAX 100

// This is where the grid is initialized.
// The main init function is at the end, and all functions are available to it.

// Uniform Sphere Model
void uniformDensity(double * Density, int Nr, int Ntheta){
  for (int i=0;i<Nr;i++){
    for (int j=0;j<Ntheta;j++){
      Density[i*Ntheta+j] = 1.e-18;
    }
  }
} 

/* 
// Uniform finite angular size Model
void Grid::initDensity(){
  double theta0 = PI/18.; // Total openning angle is 20 degree;
  for (int i=0;i<Nr;i++){
    for (int j=0;j<Ntheta;j++){
      if (thetac[j] > PI/2-theta0 && thetac[j] < PI/2+theta0)
      {
        Density[i*Ntheta+j] = 1.e-14;
      }
      else Density[i*Ntheta+j] = 1.e-18;
    }
  }
}
*/
void uniformBnuT(double * BnuT, int Nr, int Ntheta){
  for (int i=0;i<Nr;i++){
    for (int j=0;j<Ntheta;j++){
      BnuT[i*Ntheta+j] = 1;
    }
  }
}

void uniformBz(double * Bfield, int Nr, int Ntheta){
  // Uniform B field in z direction for now.
  // Note that since this is pure axis-symmetric code, By is simply Bphi.
  for (int i=0;i<Nr;i++){
    for (int j=0;j<Ntheta;j++){
      Bfield[(i*Ntheta+j)*3 + 0] = 0; // Bx
      Bfield[(i*Ntheta+j)*3 + 1] = 0; // By
      Bfield[(i*Ntheta+j)*3 + 2] = 1.; // Bz
    }
  }
}

void uniformRGrid(double * rc, double * rl, double * rr, int Nr, double& epsDS){
  double dr = RMAX*AU/Nr;
  for(int i=0;i<Nr;i++){
    rc[i] = dr*(i+0.5);
    rl[i] = dr*i;
    rr[i] = dr*(i+1);
  }
  epsDS = dr/1e5;
}

void uniformThetaGrid(double * thetac, double * thetal, double * thetar, int Ntheta){
  for(int i=0;i<Ntheta;i++){ // Max theta is PI now. Note that it is pretty bad.
    // I'll work on a modification later. 
    thetac[i] = PI/Ntheta*(i+0.5);
    thetal[i] = PI/Ntheta*(i);
    thetar[i] = PI/Ntheta*(i+1.);
  }
}

void Grid::init(){
  uniformRGrid(rc, rl, rr, Nr, epsDS);
  uniformThetaGrid(thetac, thetal, thetar, Ntheta);
  // Angular grid is not customizable at this point since it involves
  // how the angular integration is done.

  uniformDensity(Density, Nr, Ntheta);
  uniformBnuT(BnuT, Nr, Ntheta);
  uniformBz(Bfield, Nr, Ntheta);
}


