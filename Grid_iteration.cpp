#include "Grid.h"
#include <iostream>
#include <fstream>

void Grid::iteration(){
  //std::ofstream Fout;
  //Fout.open("test/positions");

  double r0, theta0; // Initial location of the cell
  double n_phi, n_theta; // Direction of the line in question
  double r, theta; // Current location of the calculation point
  double x, y, z;
  double nx, ny, nz;
  double ds, dtau;
  double I,Q,U,V;
  int ir, it, irs, its;
  double rho, bnuT;
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
	  //std::cout<<"Starting point: "<<r/AU<<","<<theta<<std::endl;
	  //std::cout<<"Angle: "<<n_theta<<","<<n_phi<<std::endl;
	  double s=0, tau=0;
	  while (this->isInDomain(x, y, z)){ // If still in the domain
	    //Fout<<x/AU<<" "<<y/AU<<" "<<z/AU<<std::endl;
	    irs = ir; its = it; // Save location information for scattering.
	    rho = this->get_density(ir,it);
	    bnuT = this->get_bnuT(ir, it);
	    //std::cout<<"ir, it (before): "<<ir<<", "<<it<<std::endl;
	    this->moveOneCell(x, y, z, nx, ny, nz, // Parameters
	    	ds, ir, it); // Things to change
	    //std::cout<<"ir, it (after) : "<<ir<<", "<<it<<std::endl;
	    //std::cout<<"distance: "<<ds/AU<<std::endl;
	    dtau =  rho * ds * kappa_ext;

	    I += bnuT * kappa_abs/kappa_ext * 
	      (exp(-tau) - exp(-(tau+dtau))); // Try new scheme.
	    //I += rho * bnuT * ds * kappa_abs * exp(-(tau+0.5*dtau)); 
	       // Thermal emission part. Non-polarized for now.

	    this->calc_Scattering(irs,its, x,y,z, nx,ny,nz, tau, dtau, // Parameters
	    	I, Q, U, V); // Things to change
	    x -= nx*ds; y -= ny*ds; z -= nz*ds; // Opposite direction.

	    //s+=ds;
	    tau+=dtau;
	    //std::cout<<"tau: "<<tau<<std::endl;
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
	  
	  //std::cout<<s/AU<<std::endl;
	  //std::cout<<I<<" "<<Q<<" "<<U<<" "<<V<<std::endl;

	  //goto endloop;
	}
      }
    }
  }
  //endloop:
  //Fout.close();
  return;
}
