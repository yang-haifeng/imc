#include "Grid.h"

Matrix Grid::calcExtinction(int ir, int it, double x, double y, double z, 
    double nx, double ny, double nz){
  Matrix S;
  for (int i=0; i<4; i++){
    for (int j=0; j<4; j++){
      if (i==j) S[4*i+j] = kappa_ext;
      else S[4*i+j] = 0.;
    }
  }
  return S;
}
