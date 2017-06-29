#include "Grid.h"
//#include <iostream>

double Grid::get_density(int ir, int it){
  return Density[ir*Ntheta+it];
}

double Grid::get_bnuT(int ir, int it){
  return BnuT[ir*Ntheta+it];
}

bool Grid::isInDomain(double x, double y, double z){
// This function return a bool value that tells if a point is in the current
// calculation domain
  double r2 = x*x+y*y+z*z;
  if(r2<r2max && r2>r2min) return true;
  else return false;
}

bool Grid::findCell(double x, double y, double z, int &ir, int &it){
  if (!this->isInDomain(x,y,z)) return false;
  double r, theta;
  r = sqrt(x*x+y*y+z*z);
  theta = acos(z/r);
  //std::cout<<"In findCell function"<<x/AU<<" "<<y/AU<<" "<<z/AU<<" "<<r/AU<<" "<<theta/PI*180<<std::endl;

  int low, high, mid;
  if (r>rl[Nr-1]) ir=Nr-1;
  else{
    low=0; high=Nr-1; 
    while(high!=low){
      mid=(low+high)/2;
      if (r>=rl[mid] && r<rr[mid]){
        ir=mid; break;
      }
      else if (r<rl[mid]) high=mid;
      else low=mid;
    }
  }

  if (theta>thetal[Ntheta-1]) it=Ntheta-1;
  else{
    low=0; high=Ntheta-1; 
    while(high!=low){
      mid=(low+high)/2;
      if (theta>=thetal[mid] && theta<thetar[mid]){
        it=mid; break;
      }
      else if (theta<thetal[mid]) high=mid;
      else low=mid;
    }
  }
  //std::cout<<ir<<" "<<it<<std::endl;
  return true;
}
