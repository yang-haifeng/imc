#include "Grid.h"

Vector Grid::Integrate(double x0, double y0, double z0, 
    double nx,double ny,double nz, bool ScaFlag){
  int ir, it;
  int irs, its;
  this->findCell(x0,y0,z0,ir,it);

  double x=x0, y=y0, z=z0;
  double rho;
  double ds, dtau, tau=0;
  Vector S,dSe,dSs; for(int i=0; i<4; i++) S[i]=0;
  while (this->isInDomain(x,y,z)){
    irs = ir; its=it;
    rho = this->get_density(ir, it);
    this->moveOneCell(x,y,z, nx,ny,nz, ds, ir,it);
    dtau = rho * ds * kappa_ext;

    dSe = this->calcEmission(irs,its, x,y,z, nx,ny,nz);

    if(ScaFlag) dSs=this->calcScattering(irs,its, x,y,z, nx,ny,nz);
    x-=nx*ds; y-=ny*ds; z-=nz*ds;

    S += (dSe+dSs) /kappa_ext * (exp(-tau) - exp(-(tau+dtau)));

    tau+=dtau;
    if(tau>10) break;
  }
  return S;
}
