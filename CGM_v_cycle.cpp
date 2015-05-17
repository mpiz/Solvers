#include "CGM_v_cycle.h"

void CGM_v_cycle::init(int *gi_s, int *gj_s, double *di_s, double *gg_s, int n_s){
	gi = gi_s;
	gj = gj_s;
	di = di_s;
	gg = gg_s;
	n = n_s;

	r = new double [n];
	z = new double [n];
	s = new double [n];

	mul_matrix = &CGM_v_cycle::mul_matrix_main;
}

void CGM_v_cycle::init_projector(int *gi_R_s, int *gj_R_s, double *gg_R_s, int n_lvl2_s) {

	gi_R = gi_R_s;
	gj_R = gj_R_s;
	gg_R = gg_R_s;
	n_lvl1 = n;
	n = n_lvl2_s;

	delete[] r;
	delete[] z;
	delete[] s;

	r = new double [n];
	z = new double [n];
	s = new double [n];

	mul_R1 = new double [n_lvl1];
	mul_R2 = new double [n_lvl1];

	mul_matrix = &CGM_v_cycle::mul_matrix_proj;
}


double CGM_v_cycle::dot_prod(double *a, double *b){
	double d_p = 0;
	for(int i = 0; i < n; i++)
		d_p += a[i]*b[i];
	return d_p;
}

void CGM_v_cycle::mul_matrix_main(double *f, double *&x){

	for(int i = 0; i < n; i++){
		double v_el = f[i];
		x[i] = di[i]*v_el;
		for(int k = gi[i], k1 = gi[i+1]; k < k1; k++){
			int j = gj[k];
			x[i] += gg[k]*f[j];
			x[j] += gg[k]*v_el;
		}		
	}
}

void CGM_v_cycle::mul_matrix_proj(double *f, double *&x){

	//R^T
	for(int i = 0; i < n_lvl1; i++)
		mul_R1[i] = 0;

	for(int i = 0; i < n; i++) {
			for(int j = gi_R[i]; j < gi_R[i+1]; j++)
				mul_R1[gj_R[j]] += gg_R[j] * f[i];
	}

	//A
	for(int i = 0; i < n_lvl1; i++){
		double v_el = mul_R1[i];
		mul_R2[i] = di[i]*v_el;
		for(int k = gi[i], k1 = gi[i+1]; k < k1; k++){
			int j = gj[k];
			mul_R2[i] += gg[k]*mul_R1[j];
			mul_R2[j] += gg[k]*v_el;
		}		
	}

	//R
	for(int i = 0; i < n; i++) {
		x[i] = 0;
		for(int j = gi_R[i]; j < gi_R[i+1]; j++) 
			x[i] += gg_R[j] * mul_R2[gj_R[j]];
	}
}


void CGM_v_cycle::solve(double *&solution, double *right_part){

	// Параметры решателя
	int max_iter = 1000;
	double eps = 1E-20;
	rp = right_part;
	x0 = new double [n];

	for(int i = 0; i < n; i++)
		x0[i] = 0;

	(this->*mul_matrix)(x0, r);

	for(int i = 0; i < n ; i++)
		z[i] = r[i] = rp[i] - r[i];


	double alpha, betta, prod_1, prod_2;
	double discr, rp_norm;

	rp_norm = sqrt(dot_prod(rp,rp));

	prod_1 = dot_prod(r, r);

	bool end = false;

	for(int iter = 0; iter < max_iter && !end; iter++){

		discr = sqrt(dot_prod(r,r));
		if(eps < discr/rp_norm){

			(this->*mul_matrix)(z, s);

			alpha = prod_1 / dot_prod(s, z);

			for(int i = 0; i < n ; i++){
				x0[i] += alpha * z[i];
				r[i] -= alpha * s[i];
			}

			prod_2 = dot_prod(r, r);

			betta = prod_2 / prod_1;

			prod_1 = prod_2;

			for(int i = 0; i < n; i++)
				z[i] = r[i] + betta*z[i];
		}
		else
			end = true;
	}

	for(int i = 0; i < n; i++)
		solution[i] = x0[i];

	delete[] x0;

}


void CGM_v_cycle::solve_iteration(int int_num, double *solution, double *right_part){

	// Параметры решателя
	int max_iter = int_num;
	double eps = 1E-20;

	x0 = solution;
	rp = right_part;

	(this->*mul_matrix)(x0, r);

	for(int i = 0; i < n ; i++)
		z[i] = r[i] = rp[i] - r[i];


	double alpha, betta, prod_1, prod_2;
	double discr, rp_norm;

	rp_norm = sqrt(dot_prod(rp,rp));

	prod_1 = dot_prod(r, r);

	bool end = false;

	for(int iter = 0; iter < max_iter && !end; iter++){

		discr = sqrt(dot_prod(r,r));
		if(eps < discr/rp_norm){

			(this->*mul_matrix)(z, s);

			alpha = prod_1 / dot_prod(s, z);

			for(int i = 0; i < n ; i++){
				x0[i] += alpha * z[i];
				r[i] -= alpha * s[i];
			}

			prod_2 = dot_prod(r, r);

			betta = prod_2 / prod_1;

			prod_1 = prod_2;

			for(int i = 0; i < n; i++)
				z[i] = r[i] + betta*z[i];
		}
		else
			end = true;
	}
	
}


void CGM_v_cycle::clear() {
	delete[] r;
	delete[] z;
	delete[] s;
}

CGM_v_cycle::~CGM_v_cycle() {
	clear();
}
