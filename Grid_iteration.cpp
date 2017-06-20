#include "Grid.h"
#include <iostream>

#ifdef _MPI_
#include <assert.h>
#endif

void Grid::iteration(bool ScaFlag){
#ifdef _MPI_
  assert(world_size>1);
  if (my_rank==0){
    if (ScaFlag)
      std::cout<<"Iteration start."<<std::endl;
    else
      std::cout<<"Start of Zeroth iteration with no scattering."<<std::endl;
  }
#else
  if (ScaFlag)
    std::cout<<"Iteration start."<<std::endl;
  else
    std::cout<<"Start of Zeroth iteration with no scattering."<<std::endl;
#endif

  double * Stokes1;
  int NMaxStokes = Nr*Ntheta*NphiI*NthetaI;

#ifdef _MPI_
  if (my_rank==0) Stokes1 = new double [ NMaxStokes*4 ];

  int N;
  //int Ntasks = 100, Nstart, Nstop, N;
  //dN = round(double(NMaxStokes) / Ntasks);
#else
  Stokes1 = new double [ NMaxStokes*4 ];
#endif

  double r0, theta0; // Initial location of the cell
  double n_phi, n_theta; // Direction of the line in question
  double x, y, z;
  double nx, ny, nz;
  Vector S;
  int ir, it;
  int Ncal=0;
#ifdef _MPI_
  typedef struct{
    int N;
    int Stemp[4];
  } Mydata;
  Mydata data;
  MPI_Status status;

  if (my_rank == 0){ // master node
    for(int i=world_size-1; i<NMaxStokes; i++){
      MPI_Recv((char *) &data, sizeof(Mydata), MPI_BYTE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
      Stokes1[data.N*4+0] = data.Stemp[0]; Stokes1[data.N*4+1] = data.Stemp[1]; Stokes1[data.N*4+2] = data.Stemp[2]; Stokes1[data.N*4+3] = data.Stemp[3];
      data.N = i;
      if (i%100 == 0){
      std::cout<<i<<" out of "<<NMaxStokes<<" done. \r";
      std::cout.flush();
      }
      MPI_Send((char *) &data, sizeof(Mydata), MPI_BYTE, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
    }
    for(int i=0; i<world_size-1; i++){
      MPI_Recv((char *) &data, sizeof(Mydata), MPI_BYTE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
      Stokes1[data.N*4+0] = data.Stemp[0]; Stokes1[data.N*4+1] = data.Stemp[1]; Stokes1[data.N*4+2] = data.Stemp[2]; Stokes1[data.N*4+3] = data.Stemp[3];
      data.N = NMaxStokes;
      MPI_Send((char *) &data, sizeof(Mydata), MPI_BYTE, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
    }
  }
  else{ // other node
    data.N = my_rank-1;
    while(true){
      if (data.N==NMaxStokes) break;
      // Calculation
      // N = i*Ntheta*NphiI*NthetaI + j*NphiI*NthetaI + k*NthetaI + l
      N = data.N;
      int l = N%NthetaI; N /= NthetaI;
      int k = N%NphiI;   N /= NphiI;
      int j = N%Ntheta;  N /= Ntheta;
      int i = N;
      r0 = rc[i]; theta0 = thetac[j]; n_phi=phiIc[k]; n_theta=thetaIc[l];

      // The following is the main body of calculation.
      x = r0*sin(theta0); y = 0; z = r0*cos(theta0); // Current location
      ir = i; it = j;
      // (nx,ny,nz) describes the direction of light in question.
      // Note that the direction of the light is the opposite of the direction
      // of integration.
      nx = sin(n_theta)*cos(n_phi); 
      ny = sin(n_theta)*sin(n_phi); 
      nz = cos(n_theta); 

      S = this->Integrate(x,y,z, nx,ny,nz, ScaFlag);

      data.Stemp[0] = S[0]; 
      data.Stemp[1] = S[1]; 
      data.Stemp[2] = S[2]; 
      data.Stemp[3] = S[3]; 

      MPI_Send((char*) &data, sizeof(Mydata), MPI_BYTE, 0, 1, MPI_COMM_WORLD);
      MPI_Recv((char*) &data, sizeof(Mydata), MPI_BYTE, 0, 1, MPI_COMM_WORLD, &status);
    }
  }
#else
  for(int i=0; i<Nr; i++){ // i is index for radius in spacial grid
    r0 = rc[i];
    for(int j=0; j<Ntheta; j++){ // j is index for theta in spacial grid
      theta0 = thetac[j];
      for(int k=0; k<NphiI; k++){ // k is index for phi in angular grid
        n_phi=phiIc[k];
	for(int l=0; l<NthetaI; l++){ // l is index for theta in angular grid
	  n_theta=thetaIc[l];

	  // The following is the main body of calculation.
	  x = r0*sin(theta0); y = 0; z = r0*cos(theta0); // Current location
          ir = i; it = j;
          // (nx,ny,nz) describes the direction of light in question.
	  // Note that the direction of the light is the opposite of the direction
	  // of integration.
	  nx = sin(n_theta)*cos(n_phi); 
	  ny = sin(n_theta)*sin(n_phi); 
	  nz = cos(n_theta); 

	  S = this->Integrate(x,y,z, nx,ny,nz, ScaFlag);

	  Stokes1[ Ncal*4 + 0] = S[0];
	  Stokes1[ Ncal*4 + 1] = S[1];
	  Stokes1[ Ncal*4 + 2] = S[2];
	  Stokes1[ Ncal*4 + 3] = S[3];
	  
          if (Ncal%1000==0){
	    std::cout<<Ncal<<" done."<<std::endl;
	  }
	  Ncal++;
	  
	}
      }
    }
  }
#endif
  // Copy new array to past iteration. 
#ifdef _MPI_
  //MPI_Allgatherv(Stokes1, Ntasks[my_rank], MPI_DOUBLE, Stokes, Ntasks, displs, MPI_DOUBLE, MPI_COMM_WORLD);
  if (my_rank==0) {
    for(int i=0;i<Nr*Ntheta*NphiI*NthetaI*4; i++) Stokes[i] = Stokes1[i];
    delete [] Stokes1;
  }
#else
  for(int i=0;i<Nr*Ntheta*NphiI*NthetaI*4; i++) Stokes[i] = Stokes1[i];
  delete [] Stokes1;
#endif

  return;
}
