#include "Dust.h"

Dust::Dust(){
  //eps = std::complex<double> (3.78, 0.04);
  eps = std::complex<double> (11.580114687445075, 0.41914961778663073);
  re = 100.;
  //s = 1.5;
  s = 1.5;

  this->calcAM();
}

Dust::Dust(std::complex<double> eps0, double re0, double s0){
  eps = eps0;
  re = re0;
  s = s0;

  this->calcAM();
}

Dust::~Dust(){
}
