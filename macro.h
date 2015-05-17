#pragma once

#include <array>
#include <complex>

using namespace std;


typedef complex<double> dcomplex;

#ifndef matrix
#define matrix(n) array<array<double,(n)>,(n)>  //���������� ������� ������� n
#endif

#ifndef cmatrix
#define cmatrix(n) array<array<dcomplex,(n)>,(n)>  //���������� ������� ������� n
#endif

#ifndef recmatrix
#define recmatrix(n,m) array<array<double,(m)>,(n)>  //���������� ������� ������� n
#endif
