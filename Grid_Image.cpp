#include "Grid.h"
#include <fstream>

// Generate a image of Npix x Npix pixels. 
// inc is the inclination angle. Line of sight is in (sin(i), 0, cos(i)) direction.
// fName is the output file name. (default: Image.out)
// Npix default is 100.
void Grid::Image(double inc, int Npix, std::string fName){
  std::ofstream Fout;
  Fout.open(fName.c_str(), std::ios::binary);

  double nx=sin(inc), ny=0, nz=cos(inc);
  double dx, dy, dz; 
  dy=rc[Nr]*2/(Npix-1);
  dx=dy*cos(inc);
  dz=-dy*sin(inc);
  double x0, y0, z0;
  y0 = -rc[Nr];
  x0 = y0*cos(inc);
  z0 = -y0*sin(inc);
  double x,y,z;
  double I,Q,U,V;
  bool status;
  for (int i=0; i<Npix; i++){
    for (int j=0; j<Npix; j++){
      x=x0+i*dx; y=y0+j*dy; z=z0+i*dz;

      getSurface(x, y, z, nx, ny, nz, status);
      if (status) { //status==true means the point is out of the domain.
        Fout<<0.<<" "<<0.<<" "<<0.<<" "<<0.<<std::endl; // write all 0 here.
	continue; // Move on to the next point.
      }
      imageInterpolate(x, y, z, nx, ny, nz, I, Q, U, V);
      Fout<<I<<" "<<Q<<" "<<U<<" "<<V<<std::endl;
    }
  }
}

// getSurface function starts from the point (x, y, z) and moves in the direction
// (nx, ny, nz), until reaches the edge of the calculation domain.
// status=true if the initial point is already out of the calculation domain.
void Grid::getSurface(double &x, double &y, double &z, double nx, double ny, 
	double nz, bool &status){
  status=true;
  return;
}

void Grid::imageInterpolate(double x, double y, double z, double nx, double ny,
	double nz, double &I, double &Q, double &U, double &V){
  return;
}
