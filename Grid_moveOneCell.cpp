#include "Grid.h"
#include <iostream>

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

  tr     = radWall(x, y, z, -nx, -ny, -nz, R1, R2);
  ttheta = thetaWall(x, y, z, -nx, -ny, -nz, theta1, theta2);
  tphi   = phiMove(x, y, z, -nx, -ny, -nz, dphiI); // phi is different in this code,
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

  ds+=1e-5*AU; // Move a little bit beyond to make sure it doesn't stay on the wall.

  return;
}

// In all following, (nx,ny,nz) is the moving direction.

// Return the distance to the closest wall. Positive if R1 is closer, negative if
// R2 is closer.
double radWall(double x, double y, double z, double nx, double ny, double nz, double R1, double R2){
  double b, c1, c2, delta1, delta2;
  b = 2.*(x*nx+y*ny+z*nz);
  c1 = (x*x+y*y+z*z) - R1*R1;
  c2 = (x*x+y*y+z*z) - R2*R2;
  delta1 = b*b - 4*c1;
  delta2 = b*b - 4*c2;

  double solutions[4], flags[4];
  int nsol=0;
  double tsol;
  if (delta1>=0){
    tsol = (-b-sqrt(delta1))/2;
    if(tsol>0){
      solutions[nsol] = tsol;
      flags[nsol]=-1;
      nsol++;
    }
    tsol = (-b+sqrt(delta1))/2;
    if(tsol>0){
      solutions[nsol] = tsol;
      flags[nsol]=-1;
      nsol++;
    }
  }
  if (delta2>=0){
    tsol = (-b-sqrt(delta2))/2;
    if(tsol>0){
      solutions[nsol] = tsol;
      flags[nsol]=1;
      nsol++;
    }
    tsol = (-b+sqrt(delta2))/2;
    if(tsol>0){
      solutions[nsol] = tsol;
      flags[nsol]=1;
      nsol++;
    }
  }
  if (nsol==0) return 1.e50; // No solution found. Return a huge number here.
  tsol=1.e50;
  for (int i=0; i<nsol; i++){
    if (solutions[i]<fabs(tsol)){
      tsol=solutions[i]*flags[i];
    }
  }

  return tsol;
}

// Return the distance to the closest wall. Positive if theta1 is closer, 
// negative if theta2 is closer.
double thetaWall(double x, double y, double z, double nx, double ny, double nz, double theta1, double theta2){
  double a1, a2, b1, b2, c1, c2;
  double tant1, tant2;
  tant1=tan(theta1); tant2=tan(theta2);
  a1 = nx*nx + ny*ny - nz*nz*tant1*tant1;
  a2 = nx*nx + ny*ny - nz*nz*tant2*tant2;
  b1 = 2*nx*x+2*ny*y-2*nz*z*tant1*tant1;
  b2 = 2*nx*x+2*ny*y-2*nz*z*tant2*tant2;
  c1 = x*x+y*y-z*z*tant1*tant1;
  c2 = x*x+y*y-z*z*tant2*tant2;
  double delta1, delta2;
  delta1 = b1*b1-4*a1*c1;
  delta2 = b2*b2-4*a2*c2;

  double solutions[4], flags[4];
  int nsol=0;
  double tsol;
  if(theta1>1.e-3) if (delta1>=0){
    tsol = (-b1-sqrt(delta1))/2/a1;
    if(tsol>0){
      solutions[nsol] = tsol;
      flags[nsol]=-1;
      nsol++;
      //std::cout<<"sol1:"<<tsol/AU<<std::endl;
    }
    tsol = (-b1+sqrt(delta1))/2/a1;
    if(tsol>0){
      solutions[nsol] = tsol;
      flags[nsol]=-1;
      nsol++;
      //std::cout<<"sol2:"<<tsol/AU<<std::endl;
    }
  }
  if(theta2>1.e-3) if (delta2>=0){
    tsol = (-b2-sqrt(delta2))/2/a2;
    if(tsol>0){
      solutions[nsol] = tsol;
      flags[nsol]=1;
      nsol++;
      //std::cout<<"sol3:"<<tsol/AU<<std::endl;
    }
    tsol = (-b2+sqrt(delta2))/2/a2;
    if(tsol>0){
      solutions[nsol] = tsol;
      flags[nsol]=1;
      nsol++;
      //std::cout<<"sol4:"<<tsol/AU<<std::endl;
    }
  }
  if (nsol==0) return 1.e50; // No solution found. Return a huge number here.
  tsol=1.e50;
  for (int i=0; i<nsol; i++){
    if (solutions[i]<fabs(tsol)){
      tsol=solutions[i]*flags[i];
    }
  }

  return tsol;
}

// In 2D, this may not seem necessary. 
double phiMove(double x, double y, double z, double nx, double ny, double nz, double dphi){
  return 1.e50;
}
