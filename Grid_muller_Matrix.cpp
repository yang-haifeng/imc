#include "Grid.h"

// Calculate Muller's matrix, or rather phase matrix (normalized to 1).
// Apply to (I, Q, U, V) and get how much radiation per optical depth (dI,dQ,dU,dV)
// (theta, phi) is incoming radiation direction
// (nx, ny, nz) is outgoing radiation direction
// (ir, it), the cell index is also passed in case we want aligned grains
void Grid::muller_Matrix(double theta, double phi, double nx, double ny, double nz,
    double I, double Q, double U, double V,
    double &dI, double &dQ, double &dU, double &dV,
    int ir, int it){
  return;
}
