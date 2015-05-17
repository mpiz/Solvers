#include "CGMComplex.h"

void CGMComplex::init(int *gi_s, int *gj_s, dcomplex *di_s, dcomplex *gg_s, dcomplex *rp_s, int n_s){
	gi = gi_s;
	gj = gj_s;
	di = di_s;
	gg = gg_s;
	rp = rp_s;
	n = n_s;

	int m = gi[n];
	r = new dcomplex [n];
	x0 = new dcomplex [n];
	z = new dcomplex [n];
	p = new dcomplex [n];
	s = new dcomplex [n];

	L_di = new dcomplex [n];
	L_gg = new dcomplex [m];

	for(int i = 0; i < n; i++){
		L_di[i] = di[i];
		x0[i] = 0;
	}

	for(int i = 0 ; i < m ; i++)
		L_gg[i] = gg[i];
}

void CGMComplex::make_LLT_decomposition(){

	dcomplex sum_d, sum_l;

	for(int k = 0; k < n ; k++){

		sum_d = 0;
		int i_s = gi[k], i_e = gi[k+1];

		for(int i = i_s; i < i_e ; i++){

			sum_l = 0; 
			int j_s = gi[gj[i]], j_e = gi[gj[i]+1];

			for(int m = i_s; m < i; m++){
				for(int j = j_s; j < j_e; j++){
					if(gj[m] == gj[j])
					{
						sum_l += L_gg[m]*L_gg[j];
						j_s++;
					}
				}
			}
			L_gg[i] = (L_gg[i] -  sum_l)/L_di[gj[i]];

			sum_d += L_gg[i]*L_gg[i];
		}
		L_di[k] = sqrt(L_di[k] - sum_d);

	}


	for(int i = 0; i < n; i++)
		L_di[i] = 1;

	for(int i = 0; i < gi[n]; i++)
		L_gg[i] = 0;
}

dcomplex CGMComplex::dot_prod(dcomplex *a, dcomplex *b){
	dcomplex d_p = 0;
	for(int i = 0; i < n; i++)
		d_p += a[i]*b[i];
	return d_p;
}

dcomplex CGMComplex::dot_prod_cj(dcomplex *a, dcomplex *b){
	dcomplex d_p = 0;
	for(int i = 0; i < n; i++)
		d_p += conj(a[i])*b[i];
	return d_p;
}

void CGMComplex::mul_matrix(dcomplex *f, dcomplex *&x){

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

void CGMComplex::solve_L(dcomplex *f, dcomplex *&x){

	for(int k = 1, k1 = 0; k <= n; k++, k1++){
		dcomplex sum = 0;

		for(int i = gi[k1]; i < gi[k]; i++)
			sum += L_gg[i]*x[gj[i]];

		x[k1] = (f[k1] - sum)/L_di[k1];
	}

}

void CGMComplex::solve_LT(dcomplex *f, dcomplex *&x){

	for(int k = n, k1 = n-1; k > 0; k--, k1--){

		x[k1] = f[k1]/L_di[k1];
		dcomplex v_el = x[k1];

		for(int i = gi[k1]; i < gi[k]; i++)
			f[gj[i]] -= L_gg[i]*v_el;
	}
}

void CGMComplex::solve_LLT(dcomplex *f, dcomplex *&x){
	solve_L(f, x);
	solve_LT(x,x);
}

void CGMComplex::solve(dcomplex *&solution){

	// Параметры решателя
	int max_iter = 1000;
	double eps = 1E-17;

	mul_matrix(x0, r);
	make_LLT_decomposition();

	for(int i = 0; i < n ; i++)
		r[i] = rp[i] - r[i];

	solve_LLT(r, z);
	for(int i = 0; i < n; i++)
		p[i] = z[i];

	dcomplex alpha, betta, prod_1, prod_2;
	double discr, rp_norm;

	rp_norm = sqrt(dot_prod_cj(rp,rp).real());

	prod_1 = dot_prod(p, r);

	bool end = false;

	for(int iter = 0; iter < max_iter && !end; iter++){

		discr = sqrt(dot_prod_cj(r,r).real());
		if(discr/rp_norm > eps){

			mul_matrix(z, s);

			alpha = prod_1 / dot_prod(s, z);

			for(int i = 0; i < n ; i++){
				x0[i] += alpha * z[i];
				r[i] -= alpha * s[i];
			}

			solve_LLT(r, p);
			prod_2 = dot_prod(p, r);

			betta = prod_2 / prod_1;

			prod_1 = prod_2;

			for(int i = 0; i < n; i++)
				z[i] = p[i] + betta*z[i];

		}
		else
			end = true;
	}
	discr = sqrt(dot_prod(r,r).real());
	printf("Residual:\t%.3e\n", discr / rp_norm);
	for(int i = 0; i < n; i++)
		solution[i] = x0[i];

	delete[] r;
	delete[] x0;
	delete[] z;
	delete[] p;
	delete[] s;

	delete[] L_di;
	delete[] L_gg;

}