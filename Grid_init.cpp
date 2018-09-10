#include "Grid.h"

#define RMAX 100

#define Z_Bfield 2
#define Tor_Bfield 1
#define Rad_Bfield 0

// This is where the grid is initialized.
// The main init function is at the end, and all functions are available to it.

// Uniform Sphere Model
void uniformDensity(double * Density, int Nr, int Ntheta){
  for (int i=0;i<Nr;i++){
    for (int j=0;j<Ntheta;j++){
      Density[i*Ntheta+j] = 1.e-18;
    }
  }
} 

/* 
// Uniform finite angular size Model
void Grid::initDensity(){
  double theta0 = PI/18.; // Total openning angle is 20 degree;
  for (int i=0;i<Nr;i++){
    for (int j=0;j<Ntheta;j++){
      if (thetac[j] > PI/2-theta0 && thetac[j] < PI/2+theta0)
      {
        Density[i*Ntheta+j] = 1.e-14;
      }
      else Density[i*Ntheta+j] = 1.e-18;
    }
  }
}
*/
void uniformBnuT(double * BnuT, int Nr, int Ntheta){
  for (int i=0;i<Nr;i++){
    for (int j=0;j<Ntheta;j++){
      BnuT[i*Ntheta+j] = 1;
    }
  }
}

void uniformBfield(double * Bfield, int Nr, int Ntheta, int Option){
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
  epsDS = dr/1e5;
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
  thetal[0] = 0.; thetar[0] = PI/2.-theta0; thetac[0] = 0.5*(thetal[0]+thetar[0]);
  thetal[Ntheta-1] = PI/2.+theta0; thetar[Ntheta-1] = PI;
  thetac[Ntheta-1] = 0.5*(thetal[Ntheta-1]+thetar[Ntheta-1]);

  double dtheta = 2*theta0/(Ntheta-2);
  for(int i=1; i<Ntheta-1; i++){
    thetal[i] = thetar[i-1];
    thetar[i] = thetal[i]+dtheta;
    thetac[i] = 0.5*(thetal[i]+thetar[i]);
  }
}

void ConicDensity(double * Density, int Nr, int Ntheta){
  for (int i=0;i<Nr;i++){
    Density[i*Ntheta+0] = 0;
    for (int j=1;j<Ntheta-1;j++){
      Density[i*Ntheta+j] = 1.e-15;
    }
    Density[i*Ntheta+Ntheta-1] = 0;
  }
}

void HLTau(double * rc, double * rl, double * rr, int Nr, double &epsDS,
	double * thetac, double * thetal, double * thetar, int Ntheta,
	double * Density, double * BnuT, double * Bfield){
  uniformRGrid(rc, rl, rr, Nr, epsDS);
  //ConicThetaGrid(0.1, thetac, thetal, thetar, Ntheta);
  uniformThetaGrid(thetac, thetal, thetar, Ntheta);

  double r_max = 200*AU; double Rc = 79*AU;
  double T0 = 70; double H0 = 32.*AU;
  double Ts0 = 400.; double rs0=3.*AU; double R0 = 10.*AU;
  double p = 1.064; double q = 0.43;
  //double rho0 = 1.e-13;
  double rho0 = 1.e-16;
  //double rho0 = 5.e-15;
  //double lambda = 0.087;
  double lambda = 0.125;

  for (int i=0; i<Nr; i++){
    Density[i*Ntheta+0] = 0.;
    BnuT[i*Ntheta+0] = 0.;

    for (int j=1; j<Ntheta-1; j++){
    //for (int j=0; j<Ntheta; j++){
      double R = rc[i]; double theta = PI/2.-thetac[j];
      double thetSH = H0*pow(R/Rc, 1.5-q/2)/R;
      double f1 = pow(R/Rc, -p);
      double f2 = exp(-pow(R/Rc, 3.5-p-q/2.));
      double f3 = exp(-theta*theta/thetSH/thetSH);
      double rho = rho0 * f1 * f2 * f3;
      //std::cout<<thetSH<<std::endl;
      //std::cout<<rho<<' '<<f1<<' '<<f2<<' '<<f3<<std::endl;
      Density[i*Ntheta+j] = rho;
      //Density[i*Ntheta+j] = 1.e-16;
      //Density[i*Ntheta+j] = rho0;
      //Density[i*Ntheta+j] = 1.e-18;
      
      double W = exp(-pow(theta/3/thetSH, 2));
      double Td = W*T0*pow(0.5*R0/R, q) + (1-W)*Ts0*pow(rs0/R, q);

      BnuT[i*Ntheta+j] = planck_bnuT(Td, con_c/lambda);
      //BnuT[i*Ntheta+j] = 1.;
    }

    Density[i*Ntheta+Ntheta-1] = 0.;
    BnuT[i*Ntheta+Ntheta-1] = 0.;
  }

  //uniformBfield(Bfield, Nr, Ntheta, Rad_Bfield);
  //uniformBfield(Bfield, Nr, Ntheta, Tor_Bfield);
  uniformBfield(Bfield, Nr, Ntheta, Z_Bfield);
}

void LP74(double * rc, double * rl, double * rr, int Nr, double &epsDS,
	double * thetac, double * thetal, double * thetar, int Ntheta,
	double * Density, double * BnuT, double * Bfield){
  uniformRGrid(rc, rl, rr, Nr, epsDS);
  //ConicThetaGrid(0.1, thetac, thetal, thetar, Ntheta);
  uniformThetaGrid(thetac, thetal, thetar, Ntheta);

  double gamma = 0.3; 

  double r_max = 200*AU; double Rc = 100*AU;
  double T0 = 70; double H0 = 12.*AU * 0.25;
  double Ts0 = 400.; double rs0=3.*AU; double R0 = 10.*AU;
  double q = 0.43; double p = gamma + 1.5 - q/2.; 

  double Sigma_c = 25; // g/cm^2
  double rho0 = Sigma_c / 100. / sqrt(2*PI) / H0;
  std::cout<<"rho0 = "<<rho0<<std::endl;
  //double rho0 = 5.e-15;
  double lambda = 0.087;
  //double lambda = 0.125;

  for (int i=0; i<Nr; i++){
    Density[i*Ntheta+0] = 0.;
    BnuT[i*Ntheta+0] = 0.;

    for (int j=1; j<Ntheta-1; j++){
    //for (int j=0; j<Ntheta; j++){
      double R = rc[i]; double theta = PI/2.-thetac[j];
      double thetSH = H0*pow(R/Rc, 1.5-q/2)/R;
      double f1 = pow(R/Rc, -p);
      double f2 = exp(-pow(R/Rc, 3.5-p-q/2.));
      double f3 = exp(-theta*theta/thetSH/thetSH);
      double rho = rho0 * f1 * f2 * f3;
      //std::cout<<thetSH<<std::endl;
      //std::cout<<rho<<' '<<f1<<' '<<f2<<' '<<f3<<std::endl;
      Density[i*Ntheta+j] = rho;
      //Density[i*Ntheta+j] = 1.e-16;
      //Density[i*Ntheta+j] = rho0;
      //Density[i*Ntheta+j] = 1.e-18;
      
      double W = exp(-pow(theta/3/thetSH, 2));
      double Td = W*T0*pow(0.5*R0/R, q) + (1-W)*Ts0*pow(rs0/R, q);

      BnuT[i*Ntheta+j] = planck_bnuT(Td, con_c/lambda);
      //BnuT[i*Ntheta+j] = 1.;
    }

    Density[i*Ntheta+Ntheta-1] = 0.;
    BnuT[i*Ntheta+Ntheta-1] = 0.;
  }

  uniformBfield(Bfield, Nr, Ntheta, Rad_Bfield);
}

void beta_Pic(double * rc, double * rl, double * rr, int Nr, double &epsDS,
	double * thetac, double * thetal, double * thetar, int Ntheta,
	double * Density, double * BnuT, double * Bfield){
  uniformRGrid(rc, rl, rr, Nr, epsDS);
  //ConicThetaGrid(0.1, thetac, thetal, thetar, Ntheta);
  uniformThetaGrid(thetac, thetal, thetar, Ntheta);
  //ConicDensity(Density, Nr, Ntheta);

  for (int i=0; i<Nr; i++){
    //for (int j=1; j<Ntheta-1; j++){
    for (int j=0; j<Ntheta; j++){
      double R = rc[i]; double theta = PI/2. - thetac[j];
      BnuT[i*Ntheta+j] = pow(R/50/AU, -0.5);

      double z = R*theta;
      if (fabs(z)<10*AU){
        Density[i*Ntheta+j] = 1.e-15;
      }
      else Density[i*Ntheta+j] = 0.;
    }
  }
  uniformBfield(Bfield, Nr, Ntheta, Tor_Bfield);
}

void IRS63(double * rc, double * rl, double * rr, int Nr, double &epsDS,
	double * thetac, double * thetal, double * thetar, int Ntheta,
	double * Density, double * BnuT, double * Bfield){
  uniformRGrid(rc, rl, rr, Nr, epsDS);
  uniformThetaGrid(thetac, thetal, thetar, Ntheta);

  double Sigma0 = 0.18387776826;
  double Rout = 100*AU;
  double T0 = 30;
  double lambda = 0.125;

  for (int i=0; i<Nr; i++){
    double R = rc[i];

    double H = 0.17 * R * pow((R/Rout), 0.17);

    for (int j=0; j<Ntheta; j++){
      double theta = PI/2.-thetac[j];
      double z=R*sin(theta);

      Density[i*Ntheta+j] = Sigma0 / (R/Rout) / (2*PI) / H * exp(-0.5*(z*z/H/H));

      BnuT[i*Ntheta+j] = planck_bnuT(T0, con_c/lambda);
      //BnuT[i*Ntheta+j] = 1.;
    }
  }

  uniformBfield(Bfield, Nr, Ntheta, Rad_Bfield);
}

void L1527(double * rc, double * rl, double * rr, int Nr, double &epsDS,
	double * thetac, double * thetal, double * thetar, int Ntheta,
	double * Density, double * BnuT, double * Bfield){
  uniformRGrid(rc, rl, rr, Nr, epsDS);
  //ConicThetaGrid(0.1, thetac, thetal, thetar, Ntheta);
  uniformThetaGrid(thetac, thetal, thetar, Ntheta);

  double Rin = 0.1*AU, Rout = 84*AU;
  double T1 = 403.5;
  double q=0.5;
  double gas_to_dust = 100.;
  double Mdisk = 1.3e-2 * 2e33; // In solar mass
  double p = 1.7; 
  double Sdamp = 0.19;
  double H1 = 0.11*AU;
  double h = 1.2;

  double lambda = 0.087;

  for (int i=0; i<Nr; i++){
    double R = rc[i];
    double SigR = (2-p)*Mdisk
                 / 2/PI/( pow(Rout, 2-p) - pow(Rin, 2-p) ) 
		 * pow(R, -p);
    if(R>Rout) SigR *= Sdamp;
    double H = H1 * pow(R/AU, h);
    double T = T1 * pow(R/AU, -q);

    for (int j=0; j<Ntheta; j++){
      double theta = PI/2.-thetac[j];
      double z=R*theta;
      //std::cout<<thetSH<<std::endl;
      //std::cout<<rho<<' '<<f1<<' '<<f2<<' '<<f3<<std::endl;
      Density[i*Ntheta+j] = SigR/(sqrt(2*PI)*H) * exp(-pow(z/H, 2)/2) / gas_to_dust;
      //Density[i*Ntheta+j] = 1.e-16;

      //BnuT[i*Ntheta+j] = planck_bnuT(T, con_c/lambda);
      BnuT[i*Ntheta+j] = 1.;
    }
  }

  uniformBfield(Bfield, Nr, Ntheta, Tor_Bfield);
}



void Grid::init(){

  //kappa_ext = 1.7424; 
  //double albedo = 0.68707;
  //double albedo = 0.;
  //kappa_sca = kappa_ext*albedo;
  //kappa_abs = kappa_ext - kappa_sca;

  // Band 7; 100 micron
  //kappa_abs = 0.709422E+00;
  //kappa_sca = 0.176029E+01;
  // Band 4; 100 micron
  //kappa_abs = 0.136446E+00;
  //kappa_sca = 0.505130E-01;
  // Band 3; 100 micron
  //kappa_abs = 0.717688E-01;
  //kappa_sca = 0.113965E-01;

  // 1250
  kappa_abs = 0.435373E+00;
  kappa_sca = 0.247879E+01;
  // 1300
  //kappa_abs = 0.351123E+00;
  //kappa_sca = 0.620750E-01;

  //kappa_ext = kappa_abs + kappa_sca;
  //
  //HLTau(rc, rl, rr, Nr, epsDS,
  //      thetac, thetal, thetar, Ntheta,
  //      Density, BnuT, Bfield);

  //beta_Pic(rc, rl, rr, Nr, epsDS,
  //      thetac, thetal, thetar, Ntheta,
  //      Density, BnuT, Bfield);

  //kappa_abs = 0.031; // Why so small? 
  //kappa_sca = 0.;
  kappa_ext = kappa_abs+kappa_sca;
  IRS63(rc, rl, rr, Nr, epsDS,
        thetac, thetal, thetar, Ntheta,
        Density, BnuT, Bfield);

  //uniformRGrid(rc, rl, rr, Nr, epsDS);
  //uniformThetaGrid(thetac, thetal, thetar, Ntheta);
  //ConicThetaGrid(PI/18., thetac, thetal, thetar, Ntheta);

  // Angular grid is not customizable at this point since it involves
  // how the angular integration is done.

  //uniformDensity(Density, Nr, Ntheta);
  //ConicDensity(Density, Nr, Ntheta);
  //uniformBnuT(BnuT, Nr, Ntheta);
  //uniformBfield(Bfield, Nr, Ntheta, Rad_Bfield);
}


