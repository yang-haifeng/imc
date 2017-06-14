#include "Grid.h"
#include <iostream>
#include <fstream>

using namespace std;

void Grid::Test(){
  double theta, phi;
  double theta0, phi0;
  double nx, ny, nz;
  double i,q,u,v;
  double di,dq,du,dv;

  int N = 10;
  //dthetaI=double(N); dphiI=double(N);
  //dthetaI=2./N; dphiI=2*PI/N;
  dthetaI=PI/N; dphiI=2*PI/N;

  for (int k1=0; k1<101; k1++){
   for (int l1=0; l1<101; l1++){
    theta0 = PI/100*k1;
    phi0 = 2*PI/100*l1;
    nx = sin(theta0) * cos(phi0);
    ny = sin(theta0) * sin(phi0);
    nz = cos(theta0);
    i=q=u=v=0.;
    for (int k=0; k<N; k++){
     for (int l=0; l<N; l++){
      theta = PI/N*(k+0.5);
      phi = 2*PI/N*(l+0.5);
    //for (int k=0; k<N+1; k++){
    // for (int l=0; l<N+1; l++){
    //  theta = PI/N*k;
    //  phi = 2*PI/N*l;
     //double mu = 2./N*(k+0.5) - 1;
     //theta = acos(mu);
      this->muller_Matrix(theta, phi, nx,ny,nz, 0., 0,0);
      //if (k==0 || k==N-1) {
      //  double factor = 0.75;
      //  di*=factor;dq*=factor;du*=factor;dv*=factor;
      //}
      //if (l==0 || l==N) {di*=0.5;dq*=0.5;du*=0.5;dv*=0.5;}
      i+=di*dthetaI*dphiI*sin(theta);
      q+=dq*dthetaI*dphiI*sin(theta);
      u+=du*dthetaI*dphiI*sin(theta);
      v+=dv*dthetaI*dphiI*sin(theta);
     }
    }
    cout<<i<<" "<<q<<" "<<u<<" "<<v<<endl;
   }
  }
}
