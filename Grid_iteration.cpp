#include "Grid.h"
#include <iostream>

void Grid::iteration(){
  double r0, theta0; // Initial location of the cell
  double n_phi, n_theta; // Direction of the line in question
  double r, theta; // Current location of the calculation point
  double x, y, z;
  double nx, ny, nz;
  double ds;
  double I,Q,U,V;
  int ir, it;
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
          // (nx,ny,nz) describes the direction of light in question.
	  // Note that the direction of the light is the opposite of the direction
	  // of integration.
	  nx = sin(n_theta)*cos(n_phi); 
	  ny = sin(n_theta)*sin(n_phi); 
	  nz = cos(n_theta); 
	  I=0; Q=0; U=0; V=0;
	  //std::cout<<"Starting point: "<<r<<","<<theta<<std::endl;
	  //std::cout<<"Angle: "<<n_theta<<","<<n_phi<<std::endl;
	  double s=0;
	  while (this->isInDomain(x, y, z)){ // If still in the domain
	    this->moveOneCell(x, y, z, nx, ny, nz, // Parameters
	    	ds, ir, it); // Things to change
	    I+=this->get_density(ir, it) * this->get_bnuT(ir, it) * ds * kappa_abs; // Thermal emission part. Non-polarized for now.
	    this->calc_Scattering(ir, it, x, y, z, nx, ny, nz, ds, // Parameters
	    	I, Q, U, V); // Things to change
	    x -= nx*ds; y -= ny*ds; z -= nz*ds; // Opposite direction.

	    s+=ds;
	  }
	  Stokes[i*Ntheta*NphiI*NthetaI*4 + j * NphiI*NthetaI*4 + k*NthetaI*4 + l*4
	  	+ 0] = I;
	  Stokes[i*Ntheta*NphiI*NthetaI*4 + j * NphiI*NthetaI*4 + k*NthetaI*4 + l*4
	  	+ 1] = Q;
	  Stokes[i*Ntheta*NphiI*NthetaI*4 + j * NphiI*NthetaI*4 + k*NthetaI*4 + l*4
	  	+ 2] = U;
	  Stokes[i*Ntheta*NphiI*NthetaI*4 + j * NphiI*NthetaI*4 + k*NthetaI*4 + l*4
	  	+ 3] = V;
	  
	  std::cout<<s/AU<<std::endl;

	  goto endloop;
	}
      }
    }
  }
  endloop:
  return;
}
