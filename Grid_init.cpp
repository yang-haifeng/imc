#include "Grid.h"
#include <iostream>

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

void Grid::init(){
  uniformDensity(Density, Nr, Ntheta);
  uniformBnuT(BnuT, Nr, Ntheta);
  uniformBz(Bfield, Nr, Ntheta);
}


