#include "CGMComplex_VC.h"
#include <cstdio>

void CGMComplex_VC::init(int *gi_s, int *gj_s, dcomplex *di_s, dcomplex *gg_s, int n_s){
	gi = gi_s;
	gj = gj_s;
	di = di_s;
	gg = gg_s;
	n = n_s;

	r = new dcomplex [n];
	w = new dcomplex [n];
	s = new dcomplex [n];
	p = new dcomplex [n];
	t = new dcomplex [n];
	u = new dcomplex [n];
}

dcomplex CGMComplex_VC::dot_prod(dcomplex* a, dcomplex* b) {
	dcomplex res;
	for(int i = 0; i < n; i++) {
		res += conj(a[i]) * b[i];
	}

	return res;
}

dcomplex CGMComplex_VC::dot_prod_nocj(dcomplex* a, dcomplex* b) {
	dcomplex res;
	for(int i = 0; i < n; i++) {
		res += a[i] * b[i];
	}

	return res;
}

void CGMComplex_VC::mul_matrix(dcomplex *f, dcomplex *&x){

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

void CGMComplex_VC::solve(dcomplex *&solution, dcomplex *rp, double gamma) {

	double eps = gamma;
	int max_iter = (int)sqrt(n);


	dcomplex alpha, beta, alpha1, alpha2; 
	double rp_norm = sqrt(dot_prod(rp, rp).real());
	x0 = new dcomplex [n];

	mul_matrix(solution, t);
	for(int i = 0; i < n; i++) {
		x0[i] = solution[i];
		r[i] = rp[i] - t[i];
		w[i] = r[i] / di[i];
		p[i] = 0;

	}

	alpha1 = dot_prod_nocj(r, w);
	beta = 0;
	
	bool not_end = true;
	double discr;
	
	int iter;
	for(iter = 0; iter < max_iter && not_end; iter++) {
		discr = sqrt(dot_prod(r,r).real());

		if(iter%50 == 0)
			printf("CCGMVC Residual:\t%5d\t%.3e\r", iter, discr / rp_norm);

		if(discr / rp_norm > eps) {
			for(int i = 0; i < n; i++)
				p[i] = w[i] + beta * p[i];
			mul_matrix(p, u);
			alpha2 = dot_prod_nocj(u, p);
			alpha = alpha1 / alpha2;

			for(int i = 0; i < n; i++) {
				x0[i] += alpha * p[i];
				r[i] -= alpha * u[i];
				w[i] = r[i] / di[i];
			}


			alpha2 = alpha1;
			alpha1 = dot_prod_nocj(r, w);
			beta = alpha1 / alpha2;

	
		}
		else
			not_end = false;
	}
	printf("CCGMVC Residual:\t%5d\t%.3e\n", iter, discr / rp_norm);
	for(int i = 0; i < n; i++)
		solution[i] = x0[i];
	delete[] x0;
}

