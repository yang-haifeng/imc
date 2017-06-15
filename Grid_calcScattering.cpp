#include "Grid.h"
#include <iostream>
#include <fstream>

// This function calculate the scattering at cell (ir, it)
// We will loop through all the informations we have about this cell. And 
// calculate the contribution of those radiations towards light propagating
// in (nx, ny, nz) direction through scattering.
// (x, y, z) is needed to determine dphi, since this is a 2D code.
Vector Grid::calcScattering(int ir, int it, double x, double y, double z,  
           double nx, double ny, double nz){
  //std::ofstream Fout;
  //Fout.open("new.dat");
  double rho = this->get_density(ir,it);
  Vector S, dS, Sin; for (int i=0;i<4;i++) S[i]=0.;
  Matrix M;
  double phiI, thetaI;
  double phi=atan2(y, x);
  for (int i=0; i<NphiI; i++){
    phiI = phiIc[i];
    for (int j=0; j<NthetaI; j++){
      thetaI = thetaIc[j];

    //Sin[0] = Stokes[ir*Ntheta*NphiI*NthetaI*4 + it*NphiI*NthetaI*4 
    //  + i*NthetaI*4 + j*4 + 0];
    //Sin[1] = Stokes[ir*Ntheta*NphiI*NthetaI*4 + it*NphiI*NthetaI*4 
    //  + i*NthetaI*4 + j*4 + 1];
    //Sin[2] = Stokes[ir*Ntheta*NphiI*NthetaI*4 + it*NphiI*NthetaI*4 
    //  + i*NthetaI*4 + j*4 + 2];
    //Sin[3] = Stokes[ir*Ntheta*NphiI*NthetaI*4 + it*NphiI*NthetaI*4 
    //  + i*NthetaI*4 + j*4 + 3];
      Sin[0] = 1.; Sin[1]=0.; Sin[2]=0.; Sin[3]=0.;
      
      M = muller_Matrix(thetaI, phiI+phi, // Incoming light direction
        nx, ny, nz, // Outgoing light direction
	phi, // We need to know the phi location of the point
	ir, it); // (ir, it) is passed for future implimentation of aligned grains

      //dS = M*Sin;
      //S += dS * sin(thetaI) * dthetaI * dphiI;
      S += M * Sin * sin(thetaI) * dthetaI * dphiI;

      //Fout<<dS[0]<<" "<<dS[1]<<" "<<dS[2]<<" "<<dS[3]<<std::endl;
    }
  }

  std::cout<<S[0]<<" "<<S[1]<<" "<<S[2]<<" "<<S[3]<<std::endl;
  //Fout.close();

  return S;
}
