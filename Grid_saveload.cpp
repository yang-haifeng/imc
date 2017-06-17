#include "Grid.h"
#include <fstream>

void Grid::saveStokes(std::string fName){
  std::ofstream Fout, Fout2;
  Fout.open(fName.c_str(), std::ios::binary);
  double temp;
  for(int i=0; i<Nr; i++){
    for(int j=0; j<Ntheta; j++){
      for(int k=0; k<NphiI; k++){
        for(int l=0; l<NthetaI; l++){
	  for (int m=0; m<4; m++){
	    temp = Stokes[i*Ntheta*NphiI*NthetaI*4 + j * NphiI*NthetaI*4 + k*NthetaI*4 + l*4 +m];
	    Fout.write( reinterpret_cast<const char*>(&temp), sizeof(double));
	  }
	}
      }
    }
  }
  Fout.close();
}

void Grid::loadStokes(std::string fName){
  std::ifstream Fin;
  Fin.open(fName.c_str(), std::ios::binary);
  double temp;
  for(int i=0; i<Nr; i++){
    for(int j=0; j<Ntheta; j++){
      for(int k=0; k<NphiI; k++){
        for(int l=0; l<NthetaI; l++){
	  for (int m=0; m<4; m++){
	  Fin.read(reinterpret_cast<char*>(&temp), sizeof(double));
	  Stokes[i*Ntheta*NphiI*NthetaI*4 + j * NphiI*NthetaI*4 + k*NthetaI*4 + l*4 +m] = temp;
	  }
	}
      }
    }
  }
  Fin.close();
}
