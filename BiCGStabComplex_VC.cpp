#include "BiCGStabComplex_VC.h"
#include <cstdio>

void BiCGStabComplex_VC::init(int *gi_s, int *gj_s, dcomplex *di_s, dcomplex *gg_s, int n_s){
	gi = gi_s;
	gj = gj_s;
	di = di_s;
	gg = gg_s;
	n = n_s;

	r = new dcomplex [n];
	v = new dcomplex [n];
	s = new dcomplex [n];
	p = new dcomplex [n];
	t = new dcomplex [n];
	r2 = new dcomplex [n];
}

dcomplex BiCGStabComplex_VC::dot_prod(dcomplex* a, dcomplex* b) {
	dcomplex res;
	for(int i = 0; i < n; i++) {
		res += conj(a[i]) * b[i];
	}

	return res;
}

dcomplex BiCGStabComplex_VC::dot_prod_nocj(dcomplex* a, dcomplex* b) {
	dcomplex res;
	for(int i = 0; i < n; i++) {
		res += a[i] * b[i];
	}

	return res;
}

void BiCGStabComplex_VC::mul_matrix(dcomplex *f, dcomplex *&x){

	for(int i = 0; i < n; i++){
		dcomplex v_el = f[i];
		x[i] = di[i]*v_el;
		for(int k = gi[i], k1 = gi[i+1]; k < k1; k++){
			int j = gj[k];
			x[i] += gg[k]*f[j];
			x[j] += gg[k]*v_el;
		}		
	}
}

void BiCGStabComplex_VC::solve(dcomplex *&solution, dcomplex *rp, double gamma) {

	double eps = gamma;
	int max_iter = (int) 100 * sqrt(n);


	dcomplex omega, alpha, ro, ro_prev, beta; 
	double rp_norm = sqrt(dot_prod(rp, rp).real());
	x0 = new dcomplex [n];

	mul_matrix(solution, t);
	for(int i = 0; i < n; i++) {
		x0[i] = solution[i];
		r[i] = rp[i] - t[i];
		r2[i] = r[i];
		v[i] = 0;
		p[i] = r[i];
	}
	ro = alpha = omega = 1;

	dcomplex a1 = 0, a2 = 0, a3;

	bool not_end = true;
	double discr = 2 * sqrt(dot_prod(r,r).real()), discr_prev;
	int iter;
	for(iter = 0; iter < max_iter && not_end; iter++) {
		discr_prev = discr;
		discr = sqrt(dot_prod(r,r).real());

		if(iter%50 == 0)
			printf("BiCGStabVC Residual:\t%5d\t%.3e\t%.3e\r", iter, discr / rp_norm, gamma);

		if(discr / rp_norm > eps) {
			mul_matrix(p, v);
			a1 = dot_prod_nocj(r, r2);
			a2 = dot_prod_nocj(v, r2);
			alpha = a1 / a2;
			for(int i = 0; i < n; i++)
				s[i] = r[i] - alpha * v[i];

			mul_matrix(s, t);
			a2 = dot_prod_nocj(s, t);
			a3 = dot_prod_nocj(t, t);
			omega = a2 / a3;

			for(int i = 0; i < n; i++) {
				x0[i] += alpha * p[i] + omega * s[i];
				r[i] = s[i] - omega * t[i];

			}

			a2 = dot_prod_nocj(r, r2);
			beta = a2/a1 * alpha / omega;

			for(int i = 0; i < n; i++)
				p[i] = r[i] + beta * (p[i] - omega * v[i]);

	
		}
		else
			not_end = false;
	}
	printf("BiCGStabVC Residual:\t%5d\t%.3e\t%.3e\n", iter, discr / rp_norm, gamma);
	for(int i = 0; i < n; i++)
		solution[i] = x0[i];
	delete[] x0;
}

