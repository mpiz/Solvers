#pragma once

#include <math.h>
#include "macro.h"

class LOSComplex {
 public:
	 void init(int *gi_s, int *gj_s, dcomplex *di_s, dcomplex *gg_s, int n_s);
	 void solve(dcomplex *&solution, dcomplex *rp_s, double eps);
 private:

	 void make_LLT_decomposition();
	 void mul_matrix(dcomplex *f, dcomplex *&x);
	 void solve_L(dcomplex *f, dcomplex *&x);
	 void solve_LT(dcomplex *f, dcomplex *&x);
	 void solve_LLT(dcomplex *f, dcomplex *&x);
	 dcomplex dot_prod(dcomplex *a, dcomplex *b);
	 dcomplex dot_prod_cj(dcomplex *a, dcomplex *b);
	 double norm_sq(dcomplex* a) 

	 int n;
	 int *gi, *gj;
	 dcomplex *di, *gg, *rp, *r, *x0, *z, *p, *s;
	 dcomplex *L_di, *L_gg;
};