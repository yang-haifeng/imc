#include "Grid.h"

Vector Grid::Integrate(double x0, double y0, double z0, 
    double nx,double ny,double nz, bool ScaFlag){
  int ir, it;
  int irs, its;
  this->findCell(x0,y0,z0,ir,it);

  double x=x0, y=y0, z=z0;
  double rho, bnuT;
  double ds, dtau, tau=0;
  double I=0, Q=0, U=0, V=0;
  double dI,dQ,dU,dV;
  while (this->isInDomain(x,y,z)){
    irs = ir; its=it;
    rho = this->get_density(ir, it);
    bnuT = this->get_bnuT(ir, it);
    this->moveOneCell(x,y,z, nx,ny,nz, ds, ir,it);
    dtau = rho * ds * kappa_ext;

    this->calcEmission(irs,its, x,y,z, nx,ny,nz, dI,dQ,dU,dV);
    I += bnuT * dI /kappa_ext * (exp(-tau) - exp(-(tau+dtau))); 
    Q += bnuT * dQ /kappa_ext * (exp(-tau) - exp(-(tau+dtau))); 
    U += bnuT * dU /kappa_ext * (exp(-tau) - exp(-(tau+dtau))); 
    V += bnuT * dV /kappa_ext * (exp(-tau) - exp(-(tau+dtau)));

    if(ScaFlag) this->calc_Scattering(irs,its, x,y,z, nx,ny,nz, tau,dtau, I,Q,U,V);
    x-=nx*ds; y-=ny*ds; z-=nz*ds;

    tau+=dtau;
    if(tau>10) break;
  }
  Vector S;
  S[0]=I; S[1]=Q; S[2]=U; S[3]=V;
  return S;
}
