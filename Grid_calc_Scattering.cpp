#include "Grid.h"

// This function calculate the scattering at cell (ir, it)
// We will loop through all the informations we have about this cell. And 
// calculate the contribution of those radiations towards light propagating
// in (nx, ny, nz) direction through scattering.
// (x, y, z) is needed to determine dphi, since this is a 2D code.
void Grid::calc_Scattering(int ir, int it, double x, double y, double z,  
           double nx, double ny, double nz, double tau, double dtau, 
	   double &I, double &Q, double &U, double &V){
  double rho = this->get_density(ir,it);
  double Iin, Qin, Uin, Vin;
  double dI=0., dQ=0., dU=0., dV=0.;
  double dIdOmega=0., dQdOmega=0., dUdOmega=0., dVdOmega=0.;
  //double nx1, ny1, nz1;
  double phiI, thetaI;
  double phi;
  phi = atan2(y, x);
  for (int i=0; i<NphiI; i++){
    phiI = phiIc[i];
    for (int j=0; j<NthetaI; j++){
      thetaI = thetaIc[j];

      Iin = Stokes[ir*Ntheta*NphiI*NthetaI*4 + it*NphiI*NthetaI*4 
        + i*NthetaI*4 + j*4 + 0];
      Qin = Stokes[ir*Ntheta*NphiI*NthetaI*4 + it*NphiI*NthetaI*4 
        + i*NthetaI*4 + j*4 + 1];
      Uin = Stokes[ir*Ntheta*NphiI*NthetaI*4 + it*NphiI*NthetaI*4 
        + i*NthetaI*4 + j*4 + 2];
      Vin = Stokes[ir*Ntheta*NphiI*NthetaI*4 + it*NphiI*NthetaI*4 
        + i*NthetaI*4 + j*4 + 3];
      
      muller_Matrix(thetaI, phiI+phi, // Incoming light direction
        nx, ny, nz, // Outgoing light direction
      	Iin, Qin, Uin, Vin, // Incoming radiation
	dIdOmega, dQdOmega, dUdOmega, dVdOmega, // Outgoing radiation
	ir, it); // (ir, it) is passed for future implimentation of aligned grains

      dI += dIdOmega * sin(thetaI) * dthetaI * dphiI;
      dQ += dQdOmega * sin(thetaI) * dthetaI * dphiI;
      dU += dUdOmega * sin(thetaI) * dthetaI * dphiI;
      dV += dVdOmega * sin(thetaI) * dthetaI * dphiI;
    }
  }

  I += dI /kappa_ext * (exp(-tau) - exp(-(tau+dtau)));
  Q += dQ /kappa_ext * (exp(-tau) - exp(-(tau+dtau)));
  U += dU /kappa_ext * (exp(-tau) - exp(-(tau+dtau)));
  V += dV /kappa_ext * (exp(-tau) - exp(-(tau+dtau)));

  return;
}
