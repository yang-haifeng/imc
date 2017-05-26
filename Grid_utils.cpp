#include "Grid.h"

double Grid::get_density(int ir, int it){
  return Density[ir*Ntheta+it];
}

double Grid::get_bnuT(int ir, int it){
  return BnuT[ir*Ntheta+it];
}

bool Grid::isInDomain(double x, double y, double z){
// This function return a bool value that tells if a point is in the current
// calculation domain
  double r2 = x*x+y*y+z*z;
  if(r2<=r2max) return true;
  else return false;
}
