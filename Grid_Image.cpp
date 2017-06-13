#include "Grid.h"
#include <iostream>
#include <fstream>

// Generate a image of Npix x Npix pixels. 
// inc is the inclination angle. Line of sight is in (sin(i), 0, cos(i)) direction.
// fName is the output file name. (default: Image.out)
// Npix default is 100.
void Grid::Image(double inc, int Npix, std::string fName){
  std::ofstream Fout;
  Fout.open(fName.c_str());

  double nx=sin(inc), ny=0, nz=cos(inc);
  double dx, dy, dz; 
  dy=rc[Nr-1]*2/(Npix-1);
  dx=dy*cos(inc);
  dz=-dy*sin(inc);
  double x0, y0, z0;
  y0 = -rc[Nr-1];
  x0 = y0*cos(inc);
  z0 = -y0*sin(inc);
  double x,y,z;
  Vector S;
  bool status;
  int ir, it;
  for (int i=0; i<Npix; i++){
    for (int j=0; j<Npix; j++){
      x=x0+i*dx; y=y0+j*dy; z=z0+i*dz;

      getSurface(x, y, z, nx, ny, nz, status, ir, it);
      if (status) { //status==true means the point is out of the domain.
        Fout<<0.<<" "<<0.<<" "<<0.<<" "<<0.<<std::endl; // write all 0 here.
	continue; // Move on to the next point.
      }

      S = this->Integrate(x,y,z, nx,ny,nz, true);

      Fout<<S[0]<<" "<<S[1]<<" "<<S[2]<<" "<<S[3]<<std::endl;
    }
  }
  Fout.close();
}

// getSurface function starts from the point (x, y, z) and moves in the direction
// (nx, ny, nz), until reaches the edge of the calculation domain.
// status=true if the initial point is already out of the calculation domain.
void Grid::getSurface(double &x, double &y, double &z, double nx, double ny, 
	double nz, bool &status, int &ir, int &it){
  if (!this->findCell(x, y, z, ir, it)){ 
    status=true; return;
  }
  int ir0, it0; double ds;
  while (this->isInDomain(x,y,z)){
    ir0=ir; it0=it;
    this->moveOneCell(x,y,z, -nx,-ny,-nz, ds, ir,it); 
    		            // moveOneCell assumes moving backward.
    x += nx*ds; y += ny*ds; z += nz*ds;
  }
  ir=ir0; it=it0;
  ds = -2e-5*AU; // move back a little bit so that (x, y, z) is in domain.
    // see moveOneCell l57 for why using this step size.
  x += nx*ds; y += ny*ds; z += nz*ds;
  status=false;
  return;
}

