#pragma once

#include <cmath>
#include "macro.h"

class BiCGStabComplex_VC {
 public:
	void init(int *gi_s, int *gj_s, dcomplex *di_s, dcomplex *gg_s, int n_s);
	void solve(dcomplex *&solution, dcomplex *rp, double gamma);

 private:
	  void mul_matrix(dcomplex *f, dcomplex *&x);
	  dcomplex dot_prod(dcomplex *a, dcomplex *b);
	  dcomplex dot_prod_nocj(dcomplex *a, dcomplex *b);

	 int n;
	 int *gi, *gj;
	 dcomplex *di, *gg;

	 dcomplex *r, *r2, *v, *s, *x0, *p, *t;

};