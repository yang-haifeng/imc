#include "Grid.h"
#include <iostream>
#include <fstream>

void progressBar(int N, int Ntot);

// Generate a image of Npix x Npix pixels. 
// inc is the inclination angle. Line of sight is in (sin(i), 0, cos(i)) direction.
// fName is the output file name. (default: Image.out)
// Npix default is 100.
void Grid::Image(double inc, int Npix, std::string fName){
#ifdef _MPI_
  if (my_rank != 0) return;
#endif
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
  int Ncount=0; 
  for (int i=0; i<Npix; i++){
    for (int j=0; j<Npix; j++){
      progressBar(Ncount, Npix*Npix); Ncount++;
      x=x0+i*dx; y=y0+j*dy; z=z0+i*dz;

      getSurface(x, y, z, nx, ny, nz, status, ir, it);
      if (status) { //status==true means the point is out of the domain.
        Fout<<0.<<" "<<0.<<" "<<0.<<" "<<0.<<std::endl; // write all 0 here.
	continue; // Move on to the next point.
      }

      //S = this->Integrate(x,y,z, nx,ny,nz, true);
      S = this->Integrate(x,y,z, nx,ny,nz, false);

      Fout<<S[0]<<" "<<S[1]<<" "<<S[2]<<" "<<S[3]<<std::endl;
    }
  }
  Fout.close();
}

Vector Grid::OnePointImage(double x0, double y0, double inc){
  double nx=sin(inc), ny=0, nz=cos(inc);
  double x, y, z;
  x = x0*cos(inc);
  y = y0;
  z = -x0*sin(inc);
  Vector S;
  bool status;
  int ir, it;
  getSurface(x, y, z, nx, ny, nz, status, ir, it);
  if (status) { //status==true means the point is out of the domain.
    S[0] = S[1] = S[2] = S[3] = 0.;
    return S; // Move on to the next point.
  }

  S = this->Integrate(x,y,z, nx,ny,nz, true);

  return S;
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
  ds = -2*epsDS; // move back a little bit so that (x, y, z) is in domain.
    // This is twice the movement in moveOneCell
  x += nx*ds; y += ny*ds; z += nz*ds;
  status=false;
  return;
}

void progressBar(int N, int Ntot){
  const int barWidth = 70;
  int pos = barWidth*N/Ntot;

  for (int i = 0; i < barWidth; ++i) {
      if (i < pos) std::cout << "=";
      else if (i == pos) std::cout << ">";
      else std::cout << " ";
  }
  std::cout << "] " << int(N*100/Ntot) << " %\r";
  std::cout.flush();
}
