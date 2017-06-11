#include "Grid.h"
#include <iostream>
using namespace std;

int main(){
  //cout<<"Hello!!"<<endl;
  Grid M;
  //cout<<M.get_density(5,5)<<endl;
  //cout<<M.get_bnuT(5,5)<<endl;
  //M.zeroIter();
  //M.saveStokes();
  std::cout<<"Generating the first image."<<std::endl;
  M.Image(PI/2.);

  //M.iteration();
  //M.saveStokes("stokes1.bin");
  //std::cout<<"Generating the second image."<<std::endl;
  //M.Image(-PI/4., 100, "image1.out");
}
