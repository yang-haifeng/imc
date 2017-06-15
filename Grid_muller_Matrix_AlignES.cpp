#include "Grid.h"
#include "utils.h"

// Note that all of the stuff here are defined with double and double *.
// I didn't bother with Eigen package in this version. 
// This might be necessary in the future when rotation of Stokes paramters
// is necessary.
double product(double e1[3], double M[3][3], double e2[3]);

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

  // Here's where electrostatic approximation comes in:
  // Sij = ei * A * ej;
  // A = diag{ a1, a1, a3} in grain frame.
  double Bx, By, Bz;
  Bx = Bfield[ (ir*Ntheta+it)*3+0 ];
  By = Bfield[ (ir*Ntheta+it)*3+1 ];
  Bz = Bfield[ (ir*Ntheta+it)*3+2 ];
  double B = sqrt(Bx*Bx+By*By+Bz*Bz);
  Bx/=B; By/=B; Bz/=B; // Normalize B field first.
  double tBx, tBy; // Rotate B field by phip, the location of the point;
  tBx = Bx*cos(phip)-By*sin(phip); tBy = Bx*sin(phip)+By*cos(phip);
  Bx = tBx; By=tBy;
  double thetaB, phiB;
  thetaB = acos(Bz); phiB = atan2(By,Bx);

  double a1 = 1.; double a3 = (1-P0)/(1+P0)*a1;
  double a[3][3];
  a[0][0] = a1+ (a3-a1) * cos(phiB)*cos(phiB)*sin(thetaB)*sin(thetaB);
  a[0][1] = (a3-a1) * sin(phiB)*cos(phiB)*sin(thetaB)*sin(thetaB);
  a[0][2] = (a3-a1) * cos(thetaB)*sin(thetaB)*cos(thetaB);
  a[1][0] = a[0][1];
  a[1][1] = a1+ (a3-a1)*sin(phiB)*sin(phiB)*sin(thetaB)*sin(thetaB);
  a[1][2] = (a3-a1) * sin(thetaB)*sin(thetaB)*cos(thetaB);
  a[2][0] = a[0][2];
  a[2][1] = a[1][2];
  a[2][2] = a1 + (a3-a1)*cos(thetaB)*cos(thetaB);

  double S11, S12, S21, S22;
  S11 = product(e1t, a, e2t); S12 = product(e1p, a, e2t);
  S21 = product(e1t, a, e2p); S22 = product(e1p, a, e2p);

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

double product(double e1[3], double M[3][3], double e2[3]){
  double s = 0;
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      s += e1[i] * M[i][j] * e2[j];
    }
  }
  return s;
}
