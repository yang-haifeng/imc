#include "Dust.h"

Dust::Dust(){
  //eps = std::complex<double> (3.78, 0.04);
  eps = std::complex<double> (11.580114687445075, 0.41914961778663073);
  re = 100.;
  //s = 1.5;
  s = 1.5;

  this->calcAM();
}

Dust::~Dust(){
}
