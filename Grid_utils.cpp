#include "Grid.h"

double Grid::get_density(int ir, int it){
  return Density[ir*Ntheta+it];
}

double Grid::get_bnuT(int ir, int it){
  return BnuT[ir*Ntheta+it];
}
