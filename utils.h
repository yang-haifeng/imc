#ifndef _UTILS_H
#define _UTILS_H

#include "math.h"
#include <array>

#define Matrix std::array <double, 16>
#define Vector std::array <double, 4>

double dot(double * a, double * b);
Matrix operator*(Matrix a, Matrix b);
Vector operator*(Matrix a, Vector b);

#endif
