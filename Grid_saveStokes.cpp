#include "Grid.h"
#include <fstream>

void Grid::saveStokes(){
  std::ofstream Fout;
  Fout.open("stokes.bin", std::ios::binary);
  for(int i=0; i<Nr; i++){
    for(int j=0; j<Ntheta; j++){
      for(int k=0; k<NphiI; k++){
        for(int l=0; l<NthetaI; l++){
	  Fout.write( (char *) &Stokes[i*Ntheta*NphiI*NthetaI*4 + j * NphiI*NthetaI*4 + k*NthetaI*4 + l*4 +0], sizeof(double));
	  Fout.write( (char *) &Stokes[i*Ntheta*NphiI*NthetaI*4 + j * NphiI*NthetaI*4 + k*NthetaI*4 + l*4 +1], sizeof(double));
	  Fout.write( (char *) &Stokes[i*Ntheta*NphiI*NthetaI*4 + j * NphiI*NthetaI*4 + k*NthetaI*4 + l*4 +2], sizeof(double));
	  Fout.write( (char *) &Stokes[i*Ntheta*NphiI*NthetaI*4 + j * NphiI*NthetaI*4 + k*NthetaI*4 + l*4 +3], sizeof(double));
	}
      }
    }
  }
}
