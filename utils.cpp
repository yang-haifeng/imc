#include "utils.h"

// Calculate the dot product of two vectors of length 3. 
// Vectors shorter than 3 will cause index error.
double dot(double * a, double * b){
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
