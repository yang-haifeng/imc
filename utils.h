#ifndef _UTILS_H
#define _UTILS_H

#include "math.h"
#include <array>
#include "typedef.h"

#define Matrix std::array <double, 16>
#define Vector std::array <double, 4>

double planck_bnuT(double T, double nu);

double dot(double * a, double * b);
Matrix operator*(Matrix a, Matrix b);
Matrix operator*(Matrix a, double b);
Matrix operator*(double b, Matrix a);
Vector operator*(Matrix a, Vector b);
Vector operator*(Vector a, double b);
Vector operator/(Vector a, double b);
Vector operator+(Vector a, Vector b);

Vector& operator+=(Vector& lhs, const Vector rhs);
Matrix& operator*=(Matrix& lhs, const double rhs);
Matrix& operator-=(Matrix& lhs, const Matrix rhs);

#endif
