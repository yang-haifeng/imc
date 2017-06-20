#include "Grid.h"
#include <iostream>

#ifdef _MPI_
#include "mpi.h"
#endif

void Grid::iteration(bool ScaFlag){
#ifdef _MPI_
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
  int Nstart, Nstop, dN; // Calculate start and stop point
  dN = round(double(NMaxStokes) / world_size);
  Nstart = dN*my_rank; Nstop = dN*(my_rank+1);
  if (my_rank == world_size-1) Nstop = NMaxStokes;

  int Ntasks[world_size]; // Prepare Ntasks for each thread
  for (int i=0; i<world_size-1; i++)
    Ntasks[i] = dN*4;
  Ntasks[world_size-1] = (NMaxStokes - dN * (world_size-1))*4;
  int displs[world_size]; // displacement
  for (int i=0; i<world_size; i++)
    displs[i] = 4*dN*i;

  Stokes1 = new double [ Ntasks[my_rank]*4 ];
#else
  Stokes1 = new double [ NMaxStokes*4 ];
#endif

  double r0, theta0; // Initial location of the cell
  double n_phi, n_theta; // Direction of the line in question
  double x, y, z;
  double nx, ny, nz;
  Vector S;
  int ir, it;
  int Ncal=0, Ncount=-1;
  for(int i=0; i<Nr; i++){ // i is index for radius in spacial grid
    r0 = rc[i];
    for(int j=0; j<Ntheta; j++){ // j is index for theta in spacial grid
      theta0 = thetac[j];
      for(int k=0; k<NphiI; k++){ // k is index for phi in angular grid
        n_phi=phiIc[k];
	for(int l=0; l<NthetaI; l++){ // l is index for theta in angular grid
	  n_theta=thetaIc[l];

#ifdef _MPI_
	  Ncount++;
          if (Ncount<Nstart || Ncount >= Nstop) continue;
#endif

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
#ifdef _MPI_
	    std::cout<<Ncal<<" done from process "<<my_rank<<std::endl;
#else
	    std::cout<<Ncal<<" done."<<std::endl;
#endif
	  }
	  Ncal++;
	  
	}
      }
    }
  }
  // Copy new array to past iteration. 
#ifdef _MPI_
  MPI_Allgatherv(Stokes1, Ntasks[my_rank], MPI_DOUBLE, Stokes, Ntasks, displs, MPI_DOUBLE, MPI_COMM_WORLD);
#else
  for(int i=0;i<Nr*Ntheta*NphiI*NthetaI*4; i++) Stokes[i] = Stokes1[i];
#endif

  delete [] Stokes1;
  return;
}
