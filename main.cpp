#include "Grid.h"
#include <iostream>
using namespace std;

#ifdef _MPI_
#include "mpi.h"
#endif

int main(){
#ifdef _MPI_
  MPI_Init(NULL, NULL);
#endif


  //cout<<"Hello!!"<<endl;
  Grid M;
  //cout<<M.get_density(5,5)<<endl;
  //cout<<M.get_bnuT(5,5)<<endl;
  M.iteration(false);
  //M.iteration(true);
  //M.saveStokes();
  //std::cout<<"Generating the first image."<<std::endl;
  M.Image(PI/2);

  //M.iteration();
  //M.saveStokes("stokes1.bin");
  //std::cout<<"Generating the second image."<<std::endl;
  //M.Image(-PI/4., 100, "image1.out");

  //cout<<"Looking at (AU, AU) point"<<endl;
  //Vector S;
  //S = M.OnePointImage(AU, AU, PI/2.);
  //cout<<S[0]<<"\t"<<S[1]<<"\t"<<S[2]<<"\t"<<S[3]<<endl;


#ifdef _MPI_
  MPI_Finalize();
#endif

}
