#include "Grid.h"

void Grid::iteration(){
  double r0, theta0; // Initial location of the cell
  double n_phi, n_theta; // Direction of the line in question
  double r, theta; // Current location of the calculation point
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
	  while (this->isInDomain(r, theta)){ // If still in the domain
	    break;
	  }
	}
      }
    }
  }
}
