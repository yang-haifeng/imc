#include "Grid.h"

#define Z_Bfield 2
#define Tor_Bfield 1
#define Rad_Bfield 0

static const double TINY_R = 1e-7*AU;
static const double TINY_THETA = 1e-7;
double RMAX;

void uniformBfield(double * Bfield, int Nr, int Ntheta, int Option);
void uniformRGrid(double * rc, double * rl, double * rr, int Nr, double& epsDS);
void logRGrid(double * rc, double * rl, double * rr, int Nr, double& epsDS);
void uniformThetaGrid(double * thetac, double * thetal, double * thetar, int Ntheta);
void ConicThetaGrid(double theta0, 
    double * thetac, double * thetal, double * thetar, int Ntheta);
double MyDensity(double r, double theta);
double MyBnuT(double, double);

void Grid::init(ParameterInput* pin){
  RMAX = pin->GetOrAddReal("model", "rmax", 100);
  uniformRGrid(rc, rl, rr, Nr, epsDS);
  //logRGrid(rc, rl, rr, Nr, epsDS);
  //uniformThetaGrid(thetac, thetal, thetar, Ntheta);
  ConicThetaGrid(PI/6., thetac, thetal, thetar, Ntheta);

  InitializeDensityWFunction(MyDensity);
  InitializeBnuTWFunction(MyBnuT);

  uniformBfield(Bfield, Nr, Ntheta, Rad_Bfield);
}

void Grid::InitializeDensityWFunction(FieldFunction_ MyDensity){
  for (int i=0;i<Nr;i++)
    for (int j=0;j<Ntheta;j++)
      Density[i*Ntheta+j] = MyDensity(rc[i], thetac[j]);
}

void Grid::InitializeBnuTWFunction(FieldFunction_ MyBnuT){
  for (int i=0;i<Nr;i++)
    for (int j=0;j<Ntheta;j++)
      BnuT[i*Ntheta+j] = MyBnuT(rc[i], thetac[j]);
}

void uniformBfield(double * Bfield, int Nr, int Ntheta, int Option){
  // Note that since this is pure axis-symmetric code, By is simply Bphi.
  for (int i=0;i<Nr;i++){
    for (int j=0;j<Ntheta;j++){
      for (int k=0; k<3; k++){
        if (k == Option)
          Bfield[(i*Ntheta+j)*3 + k] = 1.;
	else
          Bfield[(i*Ntheta+j)*3 + k] = 0;
      }
    }
  }
}

void uniformRGrid(double * rc, double * rl, double * rr, int Nr, double& epsDS){
  double dr = RMAX*AU/Nr;
  for(int i=0;i<Nr;i++){
    rc[i] = dr*(i+0.5);
    rl[i] = dr*i;
    rr[i] = dr*(i+1);
  }
  for (int i=0; i<Nr-1; i++){
    rl[i+1] -= TINY_R;
    rr[i] += TINY_R;
  }
  epsDS = dr/1e5;
}

void logRGrid(double * rc, double * rl, double * rr, int Nr, double& epsDS){
  rl[0] = 0.; rr[0] = AU + TINY_R; rc[0] = 0.5*(rl[0]+rr[0]);
  double logAU = log(AU); double logRm = log(RMAX*AU);
  double dlog = (logRm-logAU)/(Nr-1);
  for(int i=1;i<Nr;i++){
    rl[i] = rr[i-1] - 2*TINY_R;
    rr[i] = exp( dlog*(i-1) + logAU ) + TINY_R;
    rc[i] = 0.5*(rl[i]+rr[i]);
  }
  epsDS = (rr[1]-rl[1])/1e5;
}

void uniformThetaGrid(double * thetac, double * thetal, double * thetar, int Ntheta){
  for(int i=0;i<Ntheta;i++){ // Max theta is PI now. Note that it is pretty bad.
    // I'll work on a modification later. 
    thetac[i] = PI/Ntheta*(i+0.5);
    thetal[i] = PI/Ntheta*(i);
    thetar[i] = PI/Ntheta*(i+1.);
  }
}

void ConicThetaGrid(double theta0, 
    double * thetac, double * thetal, double * thetar, int Ntheta){
  thetal[0] = 0.; thetar[0] = PI/2.-theta0 + TINY_THETA; 
  thetac[0] = 0.5*(thetal[0]+thetar[0]);
  thetal[Ntheta-1] = PI/2.+theta0 - TINY_THETA; thetar[Ntheta-1] = PI;
  thetac[Ntheta-1] = 0.5*(thetal[Ntheta-1]+thetar[Ntheta-1]);

  double dtheta = 2*theta0/(Ntheta-2);
  for(int i=1; i<Ntheta-1; i++){
    thetal[i] = thetar[i-1] - 2*TINY_THETA;
    thetar[i] = thetal[i]+dtheta + 2*TINY_THETA;
    thetac[i] = 0.5*(thetal[i]+thetar[i]);
  }
}

double MyDensity(double r, double theta){
  //if (fabs(theta-PI/2)<PI/6) return 1e-15;
  //else return 0;
  double Sigma0 = 0.18387776826;
  double H = 0.17*r*pow((r/RMAX/AU),0.17);
  double z = r*cos(theta);

  double rho = Sigma0 / (r/RMAX/AU)/sqrt(2*PI)/H*exp(-0.5*(z*z/H/H));

  return rho;
}

double MyBnuT(double, double){return 1;}
