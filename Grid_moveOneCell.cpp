#include "Grid.h"

// Utility functions to find the closest wall.
double radWall(double x, double y, double z, double nx, double ny, double nz, double R1, double R2);
double thetaWall(double x, double y, double z, double nx, double ny, double nz, double theta1, double theta2);
double phiMove(double x, double y, double z, double nx, double ny, double nz, double dphi);

// This function moves the location by one cell and calculate the distance moved
void Grid::moveOneCell(double x, double y, double z, 
           double nx, double ny, double nz, 
           double &ds, int &ir, int &it){
// x,y,z is the location
// nx,ny,nz is the direction of the light path. Note that this is the opposite
// of the actual moving direction.
// ds is the returned value: distance moved
// ir, it is the cell location. It gives the boundary and also will be modified.
  double r, theta, phi;
  r = pow(x*x+y*y+z*z,0.5);
  theta = acos(z/r);
  phi = acos(x/r/sin(theta));

  double tr, tphi, ttheta;

  double R1, R2;
  R1 = rl[ir]; R2 = rr[ir]; // Inner and outer boundary of the cell.
  double theta1, theta2;
  theta1 = thetal[it]; theta2 = thetar[it]; // Inner and outer boundary of the cell.

  tr     = radWall(x, y, z, nx, ny, nz, R1, R2);
  ttheta = thetaWall(x, y, z, nx, ny, nz, theta1, theta2);
  tphi   = phiMove(x, y, z, nx, ny, nz, dphiI); // phi is different in this code,
  	// because this is a 2D code. tphi is defined by changing phi by dphiI.

  double atr, atp, att;
  atr = fabs(tr); atp = fabs(ttheta); att = fabs(tphi);

  if (atr<atp){
    if (atr<att){ // atr is smallest.
      ds  = atr; 
      ir += atr/tr;
    }
    else{ // att is smallest
      ds = att;
    }
  }
  else{
    if (atp<att){ // atp is smallest
      ds  = atp;
      it += atp/ttheta;
    }
    else{ // att is smallest
      ds = att;
    }
  }

  return;
}

double radWall(double x, double y, double z, double nx, double ny, double nz, double R1, double R2){
  return 0.;
}
double thetaWall(double x, double y, double z, double nx, double ny, double nz, double theta1, double theta2){
  return 0.;
}
double phiMove(double x, double y, double z, double nx, double ny, double nz, double dphi){
  return 0.;
}
