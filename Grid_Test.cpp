#include "Grid.h"
#include <iostream>
#include <fstream>

using namespace std;

void Grid::Test(){
  double theta, phi;
  double theta0, phi0;
  double nx, ny, nz;
  Vector s, si; Matrix m;
  si[0]=1; for (int i=1; i<4; i++) si[i]=0;
  int N = 10;
  //dthetaI=double(N); dphiI=double(N);
  //dthetaI=2./N; dphiI=2*PI/N;
  dthetaI=PI/N; dphiI=2*PI/N;

  //theta0 = PI/100*50;
  //phi0 = 2*PI/100*20;
  //theta0=0; phi0=0;
  for (int k1=0; k1<101; k1++){
   for (int l1=0; l1<101; l1++){
    theta0 = PI/100*k1;
    phi0 = 2*PI/100*l1;
    nx = sin(theta0) * cos(phi0);
    ny = sin(theta0) * sin(phi0);
    nz = cos(theta0);

  //theta = PI/100*50;
  //phi = 2*PI/100*20;
  //m = this->muller_Matrix(theta, phi, nx,ny,nz, 0., 0,0);
  //s = m*si; //*dthetaI*dphiI*sin(theta);

    for (int i=0; i<4; i++) s[i] = 0;
    for (int k=0; k<N; k++){
     for (int l=0; l<N; l++){
      theta = PI/N*(k+0.5);
      phi = 2*PI/N*(l+0.5);
      m = this->muller_Matrix(theta, phi, nx,ny,nz, 0., 0,0);
      s += m*si*dthetaI*dphiI*sin(theta);
     }
    }
    cout<<s[0]<<" "<<s[1]<<" "<<s[2]<<" "<<s[3]<<endl;
   }
  }
}
