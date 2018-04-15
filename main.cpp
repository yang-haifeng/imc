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

  //Grid M;
  //M.iteration(false);
  //M.iteration(true);
  //M.saveStokes();
  //std::cout<<"Generating the first image."<<std::endl;
  //M.Image(PI/2);

  //cout<<"Looking at (AU, AU) point"<<endl;
  //Vector S;
  //S = M.OnePointImage(AU, AU, PI/2.);
  //cout<<S[0]<<"\t"<<S[1]<<"\t"<<S[2]<<"\t"<<S[3]<<endl;

  Dust D1(std::complex<double> (1.944204975015150882e+00, 1.083091351576478956e-02), 100., 1.1);
  //Dust D2(std::complex<double> (1.943494423791821690e+00, 8.070341452893656578e-03), 100., 1.5);
  //Dust D3(std::complex<double> (1.943494423791821468e+00, 5.444254152522486762e-03), 100., 1.5);
  //Dust D4(std::complex<double> (1.943494423791821690e+00, 4.158599581203851342e-03), 100., 1.5);

#ifdef _MPI_
  MPI_Finalize();
#endif

}
