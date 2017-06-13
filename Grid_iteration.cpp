#include "Grid.h"
#include <iostream>

void Grid::iteration(bool ScaFlag){
  if (ScaFlag)
    std::cout<<"Iteration start."<<std::endl;
  else
    std::cout<<"Start of Zeroth iteration with no scattering."<<std::endl;

  double r0, theta0; // Initial location of the cell
  double n_phi, n_theta; // Direction of the line in question
  double r, theta; // Current location of the calculation point
  double x, y, z;
  double nx, ny, nz;
  Vector S;
  int ir, it;
  int Ncal=0;
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

	  S = this->Integrate(x,y,z, nx,ny,nz, ScaFlag);

	  Stokes1[i*Ntheta*NphiI*NthetaI*4 + j * NphiI*NthetaI*4 + k*NthetaI*4 + l*4
	  	+ 0] = S[0];
	  Stokes1[i*Ntheta*NphiI*NthetaI*4 + j * NphiI*NthetaI*4 + k*NthetaI*4 + l*4
	  	+ 1] = S[1];
	  Stokes1[i*Ntheta*NphiI*NthetaI*4 + j * NphiI*NthetaI*4 + k*NthetaI*4 + l*4
	  	+ 2] = S[2];
	  Stokes1[i*Ntheta*NphiI*NthetaI*4 + j * NphiI*NthetaI*4 + k*NthetaI*4 + l*4
	  	+ 3] = S[3];
	  
	  if (Ncal%1000==0){
	    std::cout<<Ncal<<" done."<<std::endl;
	  }
	  Ncal++;
	  
	}
      }
    }
  }
  // Copy new array to past iteration. 
  for(int i=0;i<Nr*Ntheta*NphiI*NthetaI*4; i++) Stokes[i] = Stokes1[i];
  return;
}
