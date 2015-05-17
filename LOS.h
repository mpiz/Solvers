#pragma once

#include <math.h>

//�������� ���, ��� ���� � �������������� �������� � ����������� ������� c LU - �������������������

class LOS{
 public:
	 void init(int* s_ig, int* s_jg, double* s_di, double* s_gl, double* s_rp, int s_n); //�������������
	 void solve(double *&solution, int &its); //��������� ������� � ���������� ��������
 private:
	 int n; //����������� ����

	 //�������� �������
	 int *ig, *jg;
	 double *gu, *gl, *di;
	 double *rp;

	 //������� ��� �����������������
	 double *Uu, *Ll, *Ld;

	 void precond(); //���������� ������ L � U
	 double dot_prod(double *a, double *b); //��������� ������������
	 void mull_A(double *f, double *&x); // x = Af
	 void solve_L(double *f, double *&x); //Lx = f, ������ ���
	 void solve_U(double *f, double *&x); //Ux = f, �������� ���
};