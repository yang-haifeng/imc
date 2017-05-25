#include "Grid.h"

double Grid::get_density(int ir, int it){
  return Density[ir*Ntheta+it];
}

double Grid::get_bnuT(int ir, int it){
  return BnuT[ir*Ntheta+it];
}

bool Grid::isInDomain(double r, double theta){
// This function return a bool value that tells if a point is in the current
// calculation domain
  if(r<=rr[Nr]) return true;
  else return false;
}
