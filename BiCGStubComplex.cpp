#include "BiCGStubComplex.h"
#include <cstdio>

void BiCGStubComplex::init(int *gi_s, int *gj_s, dcomplex *di_s, dcomplex *gg_s, dcomplex *rp_s, int n_s){
	gi = gi_s;
	gj = gj_s;
	di = di_s;
	gg = gg_s;
	rp = rp_s;
	n = n_s;
}

dcomplex BiCGStubComplex::dot_prod(dcomplex* a, dcomplex* b) {
	dcomplex res;
	for(int i = 0; i < n; i++) {
		res += conj(a[i]) * b[i];
	}

	return res;
}

dcomplex BiCGStubComplex::dot_prod_nocj(dcomplex* a, dcomplex* b) {
	dcomplex res;
	for(int i = 0; i < n; i++) {
		res += a[i] * b[i];
	}

	return res;
}

void BiCGStubComplex::mul_matrix(dcomplex *f, dcomplex *&x){

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

void BiCGStubComplex::solve(dcomplex *&solution) {

	double eps = 5e-16;
	int max_iter = 10000;

	dcomplex *r, *r2, *v, *s, *x0, *p, *t;
	dcomplex omega, alpha, ro, ro_prev, beta; 
	double rp_norm = sqrt(dot_prod(rp, rp).real());

	r = new dcomplex [n];
	v = new dcomplex [n];
	s = new dcomplex [n];
	x0 = new dcomplex [n];
	p = new dcomplex [n];
	t = new dcomplex [n];
	r2 = new dcomplex [n];

	for(int i = 0; i < n; i++) {
		x0[i] = 0;
		r[i] = rp[i]; //r = b - A(x0), x0 = 0 => r = b
		r2[i] = r[i];
		v[i] = 0;
		p[i] = 0;
	}
	ro = alpha = omega = 1;


	bool not_end = true;

	for(int iter = 0; iter < max_iter && not_end; iter++) {
		double discr = sqrt(dot_prod(r,r).real());

		if(discr / rp_norm > eps) {
			ro_prev = ro;
			ro = dot_prod(r2, r);
			beta = (ro / ro_prev) * (alpha / omega);

			for(int i = 0; i < n; i++)
				p[i] = r[i] + beta * (p[i] - omega * v[i]);

			mul_matrix(p,v);
			alpha = ro / dot_prod(r2,v);
			
			for(int i = 0; i < n; i++)
				s[i] = r[i] - alpha * v[i];

			mul_matrix(s,t);
			omega = dot_prod_nocj(t,s) / dot_prod_nocj(t,t);

			for(int i = 0; i < n; i++) {
				x0[i] = x0[i] + omega * s[i] + alpha * p[i];
				r[i] = s[i] - omega * t[i];
			}
			//обновление метода
			/*if(iter%100 == 0) {
				mul_matrix(x0, t);
				for(int i = 0; i < n; i++) {
					r[i] = rp[i] - t[i];
					r2[i] = r[i];
					v[i] = 0;
					p[i] = 0;
				}
				ro = alpha = omega = 1;

			}*/

		}
		else
			not_end = false;
	}
	double discr = sqrt(dot_prod(r,r).real());
	printf("Residual:\t%.3e\n", discr / rp_norm);
	for(int i = 0; i < n; i++)
		solution[i] = x0[i];

	delete[] r;
	delete[] v;
	delete[] s;
	delete[] x0;
	delete[] p;
	delete[] t;
	delete[] r2;

}