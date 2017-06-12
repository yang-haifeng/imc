#include "Grid.h"
#include <iostream>

// Uniform Sphere Model
void Grid::initDensity(){
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

void Grid::initBnuT(){
  for (int i=0;i<Nr;i++){
    for (int j=0;j<Ntheta;j++){
      BnuT[i*Ntheta+j] = 1;
    }
  }
}

void Grid::initBfield(){
  // Uniform B field in z direction for now.
  // Note that since this is pure axis-symmetric code, By is simply Bphi.
  for (int i=0;i<Nr;i++){
    for (int j=0;j<Ntheta;j++){
      Bfield[(i*Ntheta+j)*3 + 0] = 0; // Bx
      Bfield[(i*Ntheta+j)*3 + 1] = 1.; // By
      Bfield[(i*Ntheta+j)*3 + 2] = 0; // Bz
    }
  }
}
