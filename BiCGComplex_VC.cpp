#include "BiCGComplex_VC.h"
#include <cstdio>

void BiCGComplex_VC::init(int *gi_s, int *gj_s, dcomplex *di_s, dcomplex *gg_s, int n_s){
	gi = gi_s;
	gj = gj_s;
	di = di_s;
	gg = gg_s;
	n = n_s;

	r = new dcomplex [n];
	z = new dcomplex [n];
	s = new dcomplex [n];
	p = new dcomplex [n];
	t = new dcomplex [n];
	t1 = new dcomplex [n];

}

dcomplex BiCGComplex_VC::dot_prod(dcomplex* a, dcomplex* b) {
	dcomplex res;
	for(int i = 0; i < n; i++) {
		res += conj(a[i]) * b[i];
	}

	return res;
}

dcomplex BiCGComplex_VC::dot_prod_nocj(dcomplex* a, dcomplex* b) {
	dcomplex res;
	for(int i = 0; i < n; i++) {
		res += a[i] * b[i];
	}

	return res;
}

void BiCGComplex_VC::mul_matrix(dcomplex *f, dcomplex *&x){

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


void BiCGComplex_VC::solve(dcomplex *&solution, dcomplex *rp, double gamma) {

	double eps = gamma;
	int max_iter = (int) 100 * sqrt(n);


	dcomplex dp1, dp2, alpha, beta;
	double rp_norm = sqrt(dot_prod(rp, rp).real());
	x0 = new dcomplex [n];

	mul_matrix(solution, t);
	for(int i = 0; i < n; i++) {
		x0[i] = solution[i];
		r[i] = rp[i] - t[i];
		p[i] = z[i] = s[i] = r[i];
	}

	bool not_end = true;
	double discr;
	int iter;

	dp1 = dot_prod_nocj(p, r);

	for(iter = 0; iter < max_iter && not_end; iter++) {
		discr = sqrt(dot_prod(r,r).real());

		if(iter%50 == 0)
			printf("BiCGVC Residual:\t%5d\t%.3e\t%.3e\r", iter, discr / rp_norm, gamma);

		if(discr / rp_norm > eps) {

			mul_matrix(z,t);

			alpha = dp1 / dot_prod_nocj(s, t);

			mul_matrix(s, t1);

			for(int i = 0; i < n; i++) {
				x0[i] += alpha * z[i];
				r[i] -= alpha * t[i];
				p[i] -= alpha * t1[i];

			}

			dp2 = dot_prod_nocj(p, r);
			beta = dp2 / dp1;
			dp1 = dp2;

			for(int i = 0; i < n; i++) {
				z[i] = r[i] + beta * z[i];
				s[i] = p[i] + beta * s[i];
			}

	
		}
		else
			not_end = false;
	}
	printf("BiCGVC Residual:\t%5d\t%.3e\t%.3e\n", iter, discr / rp_norm, gamma);
	for(int i = 0; i < n; i++)
		solution[i] = x0[i];
	delete[] x0;
}

