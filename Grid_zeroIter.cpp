#include "Grid.h"
#include <iostream>
#include <fstream>

void Grid::zeroIter(){
  std::cout<<"Start of Zeroth iteration with no scattering."<<std::endl;

  double r0, theta0; // Initial location of the cell
  double n_phi, n_theta; // Direction of the line in question
  double r, theta; // Current location of the calculation point
  double x, y, z;
  double nx, ny, nz;
  double ds, dtau;
  double I,Q,U,V;
  int ir, it, irs, its;
  double rho, bnuT;
  int Ncal=0;
  double dI, dQ, dU, dV;
  for(int i=0; i<Nr; i++){ // i is index for radius in spacial grid
    r0 = rc[i];
    for(int j=0; j<Ntheta; j++){ // j is index for theta in spacial grid
      theta0 = thetac[j];
      for(int k=0; k<NphiI; k++){ // k is index for phi in angular grid
        n_phi=phiIc[k];
	for(int l=0; l<NthetaI; l++){ // l is index for theta in angular grid
	  n_theta=thetaIc[l];

	  // The following is the main body of calculation.
	  r = r0; theta = theta0; // Set initial value for current location
	  x = r*sin(theta); y = 0; z = r*cos(theta); // Current location
          ir = i; it = j;
          // (nx,ny,nz) describes the direction of light in question.
	  // Note that the direction of the light is the opposite of the direction
	  // of integration.
	  nx = sin(n_theta)*cos(n_phi); 
	  ny = sin(n_theta)*sin(n_phi); 
	  nz = cos(n_theta); 
	  I=0; Q=0; U=0; V=0;
	  double tau=0;
	  while (this->isInDomain(x, y, z)){ // If still in the domain
	    irs = ir; its = it; // Save location information for scattering.
	    rho = this->get_density(ir,it);
	    bnuT = this->get_bnuT(ir, it);
	    this->moveOneCell(x, y, z, nx, ny, nz, // Parameters
	    	ds, ir, it); // Things to change
	    dtau =  rho * ds * kappa_ext;

	    this->calcEmission(irs,its, x,y,z, nx,ny,nz, dI,dQ,dU,dV);
	    I += bnuT * dI /kappa_ext * (exp(-tau) - exp(-(tau+dtau))); 
	    Q += bnuT * dQ /kappa_ext * (exp(-tau) - exp(-(tau+dtau))); 
	    U += bnuT * dU /kappa_ext * (exp(-tau) - exp(-(tau+dtau))); 
	    V += bnuT * dV /kappa_ext * (exp(-tau) - exp(-(tau+dtau))); 

	    x -= nx*ds; y -= ny*ds; z -= nz*ds; // Opposite direction.

	    tau+=dtau;
	    if (tau>10) break;
	  }
	  Stokes[i*Ntheta*NphiI*NthetaI*4 + j * NphiI*NthetaI*4 + k*NthetaI*4 + l*4
	  	+ 0] = I;
	  Stokes[i*Ntheta*NphiI*NthetaI*4 + j * NphiI*NthetaI*4 + k*NthetaI*4 + l*4
	  	+ 1] = Q;
	  Stokes[i*Ntheta*NphiI*NthetaI*4 + j * NphiI*NthetaI*4 + k*NthetaI*4 + l*4
	  	+ 2] = U;
	  Stokes[i*Ntheta*NphiI*NthetaI*4 + j * NphiI*NthetaI*4 + k*NthetaI*4 + l*4
	  	+ 3] = V;
	  
	  if (Ncal%1000==0){
	    std::cout<<Ncal<<" done."<<std::endl;
	  }
	  Ncal++;
	  
	}
      }
    }
  }
  return;
}
