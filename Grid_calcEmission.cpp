#include "Grid.h"

// In 1,2 frame:
// p = P0 sin(i)^2 / (1+P0 cos(i)^2)
// I = kappa_abs * (1+P0 cos(i)^2) / (1+P0)
// Q = kappa_abs * P0 sin(i)^2 / (1+P0)
Vector Grid::calcEmission(int ir, int it, double x, double y, double z, 
    double nx, double ny, double nz){
  double Bx, By, Bz;
  Bx = Bfield[(ir*Ntheta+it)*3+0];
  By = Bfield[(ir*Ntheta+it)*3+1];
  Bz = Bfield[(ir*Ntheta+it)*3+2];

  double B = sqrt(Bx*Bx+By*By+Bz*Bz);
  Bx/=B; By/=B; Bz/=B; // Normalize B field first.

  double phi0 = atan2(y,x); // Rotate B field by phi0, the location of the point;
  double tBx, tBy;
  tBx = Bx*cos(phi0)-By*sin(phi0); tBy = Bx*sin(phi0)+By*cos(phi0);
  Bx = tBx; By=tBy;

  double cosinc; // inclination angle with respect to B field
  cosinc = fabs(Bx*nx+By*ny+Bz*nz); 

  double theta, phi;
  theta = acos(nz); phi = atan2(ny, nx);
  double et[3]; double ep[3];
  et[0] = cos(theta)*cos(phi); et[1]=cos(theta)*sin(phi); et[2]=-sin(theta);
  ep[0] = -sin(phi);           ep[1]=cos(phi);            ep[2]=0;
  double e1[3]; //double e2[3];
  e1[0] = By*nz-Bz*ny;
  e1[1] = Bz*nx-Bx*nz;
  e1[2] = Bx*ny-By*nx;

  double cosga = e1[0]*et[0]+e1[1]*et[1]+e1[2]*et[2];
  double singa = e1[0]*ep[0]+e1[1]*ep[1]+e1[2]*ep[2];
  Vector S; double bnuT = this->get_bnuT(ir, it);
  S[0] = bnuT*kappa_abs*(1 + P0*cosinc*cosinc);// /(1+P0)
  S[1] = bnuT*kappa_abs*P0*(1-cosinc*cosinc)*(cosga*cosga - singa*singa);// /(1+P0)
  S[2] = bnuT*kappa_abs*P0*(1-cosinc*cosinc)*2*cosga*singa;// /(1+P0)
  S[3] = 0.;

  return S;
}
