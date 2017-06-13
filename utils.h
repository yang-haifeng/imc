#ifndef _UTILS_H
#define _UTILS_H

#include "math.h"

#define Matrix array <double, 16>
#define Vector array <double, 4>

using namespace std;

double dot(double * a, double * b);
Matrix operator*(Matrix a, Matrix b);
Vector operator*(Matrix a, Vector b);

#endif
