#include "LOSComplex.h"

#include "stdio.h"

void LOSComplex::init(int *gi_s, int *gj_s, dcomplex *di_s, dcomplex *gg_s, int n_s){
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
	

}


void LOSComplex::make_LLT_decomposition(){

	dcomplex sum_d, sum_l;
	int m = gi[n];
	
	L_di = new double [n];
	L_gg = new double [m];

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

dcomplex LOSComplex::dot_prod(dcomplex* a, dcomplex* b) {
	dcomplex res;
	for(int i = 0; i < n; i++) {
		res += a[i] * b[i];
	}

	return res;
}

double LOSComplex::norm_sq(dcomplex* a) {
	dcomplex res;
	for(int i = 0; i < n; i++) {
		res += abs(a[i]);
	}

	return res;
}


void LOSComplex::mul_matrix(dcomplex *f, dcomplex *&x){

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

void LOSComplex::solve_L(dcomplex *f, dcomplex *&x){

	for(int k = 1, k1 = 0; k <= n; k++, k1++){
		dcomplex sum = 0;

		for(int i = gi[k1]; i < gi[k]; i++)
			sum += L_gg[i]*x[gj[i]];

		x[k1] = (f[k1] - sum)/L_di[k1];
	}

}}
}

void LOSComplex::solve_LT(dcomplex *f, dcomplex *&x){

	for(int k = n, k1 = n-1; k > 0; k--, k1--){

		x[k1] = f[k1]/L_di[k1];
		dcomplex v_el = x[k1];

		for(int i = gi[k1]; i < gi[k]; i++)
			f[gj[i]] -= L_gg[i]*v_el;
	}
}

void LOSComplex::solve_LLT(dcomplex *f, dcomplex *&x){
	solve_L(f, x);
	solve_LT(x,x);
}

void LOS::solve(dcomplex *&solution, dcomplex *rp_s, double eps){
	
	make_LLT_decomposition();

	//Параметры решателя
	int max_iter = 100000;
	bool end_cycle = false;

	//Норма правой части, для выхода
	rp = rp_s;
	double rp_norm = norm(rp);

	//Начинаем решение
	double* x0 = new dcomplex [n]; //Приближение
	for(int i = 0; i < n; i++)
		x0[i] = 0;

	solution = new dcomplex [n];

	//r0 = L^(-1) * (f - Ax0)
	mul_matrix(x0, s);
	for(int i = 0; i < n; i++)
		s[i] = rp[i] - s[i];
	solve_L(s, r);

	//z0 = U^(-1)r0
	solve_LT(r, z);

	//p0 = L^(-1)Az0
	mul_matrix(z, s);
	solve_L(s, p);

	int iter;

	for(iter = 0; iter < max_iter && !end_cycle; iter++){
		double discr = norm(r); // Абсолютная невязка
		printf("%.3e\n", discr / rp_norm);
		if( discr / rp_norm > eps){ //Проверка условия выхода
			dcomplex dot1 = dot_prod(p, p); //(p[k-1], p[k-1])
			dcomplex alpha = dot_prod(p ,r) / dot1; //a = (p[k-1], r[k-1]) / (p[k-1], p[k-1])

			for(int i = 0; i < n; i++){
				x0[i] = x0[i] + alpha*z[i]; //x[k] = x[k-1] + a*z[k-1]
				r[i] = r[i] - alpha*p[i]; //r[k] = r[k-1] - a*p[k-1]
			}
			//betta = -(p[k-1], L^(-1)*A*U^(-1)r[k]) / (p[k-1], p[k-1])

			solve_LT(r, s); // s = U^(-1)r[k]
			mul_matrix(s, t);
			solve_L(t, t);
			dcomplex betta = - dot_prod(p, t) / dot1;

			for(int i = 0; i < n; i++){
				z[i] = s[i] + betta * z[i]; // z[k] = U^(-1)r[k] + b*z[k-1]
				p[i] = t[i] + betta * p[i]; // p[k] = L^(-1)*A*U^(-1)r[k] + b*p[k-1]
			}

			if(iter % n == 0){ //Обновление метода
				//r0 = L^(-1) * (f - Ax0)
				mul_matrix(x0, s);
				for(int i = 0; i < n; i++)
					s[i] = rp[i] - s[i];
				solve_L(s, r);

				//z0 = U^(-1)r0
				solve_LT(r, z);

				//p0 = L^(-1)Az0
				mul_matrix(z, s);
				solve_L(s, p);
			}


		}
		else{
			end_cycle = true;
		}
	}

	//solve_LT(x0, solution);
	for(int i = 0 ; i < n; i++)
		solution[i] = x0[i];

	//И отчишаем память
	delete[] x0;
	delete[] p;
	delete[] r;
	delete[] z;
	delete[] s;
	delete[] t;

	delete[] L_di;
	delete[] L_gg;

}

