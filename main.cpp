#include "Grid.h"
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
  M.iteration(false);
  //M.iteration(true);
  //M.saveStokes();
  M.Image(PI/4, 100);

  //M.loadStokes();
  //Vector S = M.OnePointImage_wID(PI/4, 4949);
  //Vector S = M.OnePointImage_wID(PI/4, 4900);
  //cout<<S[0]<<"\t"<<S[1]<<"\t"<<S[2]<<"\t"<<S[3]<<endl;

/*
double x=3.73995e+12;
double y=0;
double z=6.47778e+12;
double nx=-0.866388;
double ny=0.441447;
double nz=0.233445;
  Vector S = M.Integrate(x,y,z, nx,ny,nz, false);
  cout<<S[0]<<"\t"<<S[1]<<"\t"<<S[2]<<"\t"<<S[3]<<endl;
  */

  //M.iteration();
  //M.saveStokes("stokes1.bin");
  //std::cout<<"Generating the second image."<<std::endl;
  //M.Image(-PI/4., 100, "image1.out");

  //cout<<"Looking at (AU, AU) point"<<endl;
  //Vector S;
  //S = M.OnePointImage(AU, AU, PI/2.);
  //cout<<S[0]<<"\t"<<S[1]<<"\t"<<S[2]<<"\t"<<S[3]<<endl;

  //double x = 3.644286e13, y = -1.666297e13, z = -2.333908e13;
  //double nx = -0.8663879181530024, ny = 0.4414467431443521, nz = 0.2334453860022741;
  //double x = 36442859508008.36, y = -16662966572724.23, z = -2333908386910.082;
  //double nx = -0.8663879181530024, ny = 0.4414467431443521, nz = 0.2334453860022741;
  //int ir = 0, it = 8;
  //double ds;
  //M.moveOneCell(x,y,z,nx,ny,nz,ds,ir,it);
  //cout<<ds<<"\t"<<ir<<"\t"<<it<<endl;

  //int i, j;
  //M.findCell(3388507418593.066, -55663026433.36113, 1956619691029.254, i, j);

#ifdef _MPI_
  MPI_Finalize();
#endif

}
