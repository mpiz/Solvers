#pragma once

#include <math.h>

class CGM_v_cycle;

typedef void (CGM_v_cycle::*mm_f)(double *, double *&);

class CGM_v_cycle{
 public:
	 void init(int *gi_s, int *gj_s, double *di_s, double *gg_s, int n_s);
	 void init_projector(int *gi_R_s, int *gj_R_s, double *gg_R_s, int n_lvl2_s);
	 void solve_iteration(int int_num, double *solution, double *right_part);
	 void solve(double *&solution, double *right_part);

	 void clear();
	 ~CGM_v_cycle();

 private:

	 void mul_matrix_main(double *f, double *&x); //A*x
	 void mul_matrix_proj(double *f, double *&x); //R * A* R^T * x

	 double dot_prod(double *a, double *b);

	 mm_f mul_matrix;

	 int n;
	 int *gi, *gj;
	 double *di, *gg, *rp, *r, *x0, *z, *s;

	 int n_lvl1;
	 int *gi_R, *gj_R;
	 double *gg_R;

	 double *mul_R1, *mul_R2;
};





