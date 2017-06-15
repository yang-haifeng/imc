#include "Grid.h"
#include "utils.h"

// Note that all of the stuff here are defined with double and double *.
// I didn't bother with Eigen package in this version. 
// This might be necessary in the future when rotation of Stokes paramters
// is necessary.

// Calculate Muller's matrix, or rather phase matrix (normalized to 1).
// Apply to (I, Q, U, V) and get how much radiation per optical depth (dI,dQ,dU,dV)
// (theta, phi) is incoming radiation direction
// (nx, ny, nz) is outgoing radiation direction
// (ir, it), the cell index is also passed in case we want aligned grains
// Return the following:
//   int_Omega M S_in dOmega
Matrix Grid::muller_Matrix(double theta, double phi, 
    double nx, double ny, double nz, double phip, int ir, int it){
// This is the version with Rayleigh limit for spherical dust grains.
  double e1t[3], e1p[3], e2t[3], e2p[3]; 
  // e1t, e1p defines the two polarization axes for incoming radiation
  // e2t, e2p defines the two polarization axes for outgoing radiation
  // The Stokes parameters are defined as in Mishchenko et al. 2000: 
  //     Q = Et * Et' - Ep * Ep'

  // For (theta, phi):
  //   et = ( cos(theta) cos(phi), cos(theta) sin(phi), -sin(theta) )
  //   ep = ( -sin(phi), cos(phi), 0)
  e1t[0]=cos(theta)*cos(phi); e1t[1]=cos(theta)*sin(phi); e1t[2]=-sin(theta);
  e1p[0]= -sin(phi);          e1p[1]=cos(phi);            e1p[2]=0;

  double theta2, phi2;
  theta2 = acos(nz);
  phi2 = atan2(ny, nx);
  e2t[0]=cos(theta2)*cos(phi2); e2t[1]=cos(theta2)*sin(phi2); e2t[2]=-sin(theta2);
  e2p[0]= -sin(phi2);           e2p[1]=cos(phi2);             e2p[2]=0;

  double S11, S12, S21, S22;
  S11 = dot(e1t, e2t); S12 = dot(e1p, e2t);
  S21 = dot(e1t, e2p); S22 = dot(e1p, e2p);

  Matrix M;

  // There seems to be some difference in the definition of Stokes parameters,
  // especially Stokes U. I've changed some '-' signs and marked the places changed
  // Ref /Users/haifengyang/working/formal_sol/models.cpp l53-l68.
  M[0*4+0] = 0.5*(S11*S11 + S12*S12 + S21*S21 + S22*S22);
  M[0*4+1] = 0.5*(S11*S11 - S12*S12 + S21*S21 - S22*S22);
  M[0*4+2] = (S11*S12 + S22*S21); // Here
  M[0*4+3] = 0.;
  M[1*4+0] = 0.5*(S11*S11 + S12*S12 - S21*S21 - S22*S22);
  M[1*4+1] = 0.5*(S11*S11 - S12*S12 - S21*S21 + S22*S22);
  M[1*4+2] = (S11*S12 - S22*S21); // Here
  M[1*4+3] = 0.;
  M[2*4+0] = (S11*S21 + S22*S12); // Here
  M[2*4+1] = (S11*S21 - S22*S12); // Here
  M[2*4+2] = (S11*S22 + S12*S21);
  M[2*4+3] = 0.;
  M[3*4+0] = 0.;
  M[3*4+1] = 0.;
  M[3*4+2] = 0.;
  M[3*4+3] = (S22*S11 - S12*S21);

  M *= 3./8.*PI*kappa_sca;

  return M;
}
