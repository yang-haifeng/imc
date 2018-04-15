#include "Dust.h"
#include <iostream>
#include <assert.h>

void Dust::calcAM(){
  assert(s>1);
  double L1, L3;
  double e2,ge;
  e2 = 1-1./(s*s);
  ge = sqrt((1-e2)/e2);
  L1 = ge/2./e2*(PI/2-atan(ge)) -ge*ge/2;
  L3 = 1-2*L1;

  a1 = 4*PI*pow(re,3)*(eps-1.)/(3.+3*L1*(eps-1.));
  a3 = 4*PI*pow(re,3)*(eps-1.)/(3.+3*L3*(eps-1.));

  p0 = (a1-a3).imag()/(a1+a3).imag();

  std::cout<<p0<<std::endl;
}
