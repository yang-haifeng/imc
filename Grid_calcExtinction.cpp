#include "Grid.h"

Matrix Grid::calcExtinction(int ir, int it, double x, double y, double z, 
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

  double Cep = kappa_ext + kappa_ext*P0*cosinc*cosinc;
  double Cpp = kappa_ext * P0 * (1-cosinc*cosinc);
  double Ccpp = 0;

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
  double cos2ga = cosga*cosga - singa*singa;
  double sin2ga = 2*cosga*singa;

  Matrix M;
  M[4*0+0] = M[4*1+1] = M[4*2+2] = M[4*3+3] = Cep;
  M[4*0+1] = M[4*1+0] = Cpp*cos2ga;
  M[4*0+2] = M[4*2+0] = Cpp*sin2ga;
  M[4*0+3] = M[4*3+0] = 0.;
  M[4*1+2] = M[4*2+1] = 0.;
  M[4*1+3] = -Cpp*sin2ga; M[4*3+1] =   Cpp*sin2ga;
  M[4*2+3] = Ccpp*cos2ga; M[4*3+2] = -Ccpp*cos2ga;

/*
  Mout << Cep, -Cpp*cos2phi, -Cpp*sin2phi, 0,
-Cpp*cos2phi, Cep, 0, Ccpp*sin2phi,
-Cpp*sin2phi, 0, Cep, -Ccpp*cos2phi,
0, -Ccpp*sin2phi, Ccpp*cos2phi, Cep;
*/

  /* This is spherical dust grain extinction.
  Matrix M;
  for (int i=0; i<4; i++){
    for (int j=0; j<4; j++){
      if (i==j) M[4*i+j] = kappa_ext;
      else M[4*i+j] = 0.;
    }
  }
  */ 
  return M;
}
