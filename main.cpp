#include "Grid.h"
#include "Dust.h"
#include <iostream>
using namespace std;

#ifdef _MPI_
int world_size;
int my_rank;
#endif

int main(){
#ifdef _MPI_
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

  Grid M;
  //M.iteration(false);
  //M.iteration(true);
  //M.saveStokes();
  //M.loadStokes();
  //std::cout<<"Generating the first image."<<std::endl;
  M.Image(PI/2, 400);
  //M.Image(45./180.*PI, 100);
  //M.Image(87./180.*PI, 400);

  //cout<<"Looking at (AU, AU) point"<<endl;
  //Vector S;
  //S = M.OnePointImage(AU, AU, PI/2.);
  //cout<<S[0]<<"\t"<<S[1]<<"\t"<<S[2]<<"\t"<<S[3]<<endl;

  //Dust D;

#ifdef _MPI_
  MPI_Finalize();
#endif

}
