#include "Grid.h"
#include <iostream>
#include <iomanip>
#include <assert.h>

inline void initT(Matrix& T);

Vector Grid::Integrate(double x0, double y0, double z0, 
    double nx,double ny,double nz, bool ScaFlag){
  int ir, it;
  int irs, its;
  this->findCell(x0,y0,z0,ir,it);

  double x=x0, y=y0, z=z0;
  double rho;
  double ds, dtau, tau=0;
  Vector S,dSe,dSs; for(int i=0; i<4; i++) S[i]=0;
  Matrix M, T;
  initT(T);
  while (this->isInDomain(x,y,z)){
    irs = ir; its=it;
    rho = this->get_density(ir, it);
    this->moveOneCell(x,y,z, nx,ny,nz, ds, ir,it);
    dtau = rho * ds * kappa_ext;

    dSe = this->calcEmission(irs,its, x,y,z, nx,ny,nz);

    if(ScaFlag) dSs=this->calcScattering(irs,its, x,y,z, nx,ny,nz);
    else dSs[0]=dSs[1]=dSs[2]=dSs[3]=0.;

    M = this->calcExtinction(irs,its, x,y,z, nx,ny,nz);

    S += T * (dSe+dSs) * rho * ds; // Might be M^-1 T M instead of just T.

    T -= T * M * rho * ds;
    x-=nx*ds; y-=ny*ds; z-=nz*ds;

    //S += (dSe+dSs) /kappa_ext * (exp(-tau) - exp(-(tau+dtau)));
    //std::cout<<std::setprecision(16)<<x<<"\t"<<y<<"\t"<<z<<"\t"<<std::endl;
    //std::cout<<std::setprecision(16)<<nx<<"\t"<<ny<<"\t"<<nz<<"\t"<<ds<<std::endl;
    //std::cout<<ir<<"\t"<<it<<"\t"<<std::endl;
    //int tir, tit;
    //this->findCell(x, y, z, tir, tit);
    //if( tir!=ir || tit!=it ) {std::cout<<"Assertion failed!"<<std::endl;}

    //std::cout<<x/AU<<'\t'<<y/AU<<'\t'<<z/AU<<std::endl;
    //std::cout<<S[0]<<'\t'<<S[1]<<'\t'<<S[2]<<'\t'<<S[3]<<std::endl;
    //std::cout<<x/AU<<'\t'<<y/AU<<'\t'<<z/AU<<'\t'<<S[0]<<std::endl;

    //if (ds<1e-3*AU || ds>40*AU ){ 
    //  std::cout<<"ds wrong! "<<ds<<std::endl; 
    //  std::cout<<ir<<tir<<'\t'<<it<<tit<<std::endl;
    //  break;
    //}

    tau+=dtau;
    if(tau>10) break;
  }
  return S;
}

void initT(Matrix& T){
  for(int i=0; i<4; i++){
   for(int j=0; j<4; j++){
     if(i==j) T[i*4+j] = 1;
     else T[i*4+j] = 0;
   }
  }
}
