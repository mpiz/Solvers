#pragma once

#include <cmath>
#include "macro.h"
//в макро объ€влен тип dcomplex = std::complex<double>

class BiCGStubComplex {
 public:
	void init(int *gi_s, int *gj_s, dcomplex *di_s, dcomplex *gg_s, dcomplex *rp_s, int n_s);
	void solve(dcomplex *&solution);

 private:
	  void mul_matrix(dcomplex *f, dcomplex *&x);
	  dcomplex dot_prod(dcomplex *a, dcomplex *b);
	  dcomplex dot_prod_nocj(dcomplex *a, dcomplex *b);

	 int n;
	 int *gi, *gj;
	 dcomplex *di, *gg, *rp;



};