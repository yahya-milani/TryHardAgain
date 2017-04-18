
//#include <stdafx.h>
#include <stdio.h>
//#include <tchar.h>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

//Global Variables
const int ele_num = 96;
const int node_num = 35;
const int bound_num = 4;
const int boundary[4] = { 0, 3, 6, 7 }; //{ 2,3,5,6 };
const int DOF = 3;
const int length = DOF*(node_num - bound_num);
const double Ro = 1000;
const double m = 782.7586;
const double la = 7866.7;
const double tF = 100;
const double dt = 0.01;
const double dt_save = 1;

// Global Arrays
const double C[6][6] = { { la + 2 * m, la, la, 0, 0, 0 }, { la, la + 2 * m, la, 0, 0, 0 }, { la, la, la + 2 * m, 0, 0, 0 }, { 0, 0, 0, m, 0, 0 }, { 0, 0, 0, 0, m, 0 }, { 0, 0, 0, 0, 0, m } };
static double B_ele[ele_num][6][12];
static double V[ele_num];
static double K[length][length];
static double M[length][length];



//
void Multi(double **x, int Ix, int Jx, double **y, int Iy, int Jy);
void SpVec(double *A_sp, int *RA_sp, int *CA_sp, double *b_vec, int length_vec);
double VecVec(double *a, double *b, int n);
//double *Solve(double *A_sp, int *RA_sp, int *CA_sp, double *B, int length);
//double **Material()

int main() {
	//Nodes
	double node[4][2];
	node[0][0]=0;
	node[0][1]=0;
	node[0][0]=1;
	node[0][1]=0;
	node[0][0]=0;
	node[0][1]=1;
	node[0][0]=1;
	node[0][1]=1;

	//Create Nodes DOF ***********************************************/
	static int ID_node[4][2];
	int counter = 0;
	bool bound;
	for (int i = 0; i < 4; i++) {
		bound = 1;
		if (i==boundary[counter]) {
			counter++;
			for (int k = 0; k < 2; k++) {
				ID_node[i][k] = -1;
			}
		}
		else {
			for (int k = 0; k < 2; k++){
				ID_node[i][k] = 2*(i-counter)+ k;
			}
		}
	}

	//Read Elements***********************************************/
	int ele[2][3];
	ele[0][0] = 0;
	ele[0][1] = 1;
	ele[0][2] = 3;
	ele[1][0] = 0;
	ele[1][1] = 3;
	ele[1][2] = 2;

	//Create Elemental DOF ***********************************************/
	static int ID_ele[2][3*2];
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 2; k++) {
				ID_ele[i][j*2 + k] = ID_node[ele[i][j]][k];
			}
		}
	}

	//************************************************/
	double x[3], y[3];
	static double B[2][3][6];

	//set B matrix and voulme for every element
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 3; j++) {
			x[j] = node[ele[i][j]][0];
			y[j] = node[ele[i][j]][1];
		}

		V[i] = ((x[0]-x[2])*(y[1]-y[2])-(y[0]-y[2])*(x[1]-x[2]))/2;
		if (V[i] < 0){
			double temp = x[1];
			x[1] = x[2];
			x[2] = temp;
			temp = y[1];
			y[1] = y[2];
			y[2] = temp;
			V[i]=((x[0]-x[2])*(y[1]-y[2])-(y[0]-y[2])*(x[1]-x[2]))/2;
			for (int j =2; j<4; j++) {
				int temp2 = ID_ele[i][j];
				ID_ele[i][j] = ID_ele[i][j + 2];
				ID_ele[i][j + 2] = temp2;
			}
		}

		B[i][0][0] = (y[1]-y[2])/(2*V[i]);
		B[i][0][1] = 0;
		B[i][0][2] = (y[2]-y[0])/(2*V[i]);
		B[i][0][3] = 0;
		B[i][0][4] = (y[0]-y[1])/(2*V[i]);
		B[i][0][5] = 0;
		B[i][1][0] = 0;
		B[i][1][1] = (x[2]-x[1])/(2*V[i]);
		B[i][1][2] = 0;
		B[i][1][3] = (x[0]-x[2])/(2*V[i]);
		B[i][1][4] = 0;
		B[i][1][5] = (x[1]-x[0])/(2*V[i]);
		B[i][2][0] = (x[2]-x[1])/(2*V[i]);
		B[i][2][1] = (y[1]-y[2])/(2*V[i]);
		B[i][2][2] = (x[0]-x[2])/(2*V[i]);
		B[i][2][3] = (y[2]-y[0])/(2*V[i]);
		B[i][2][4] = (x[3]-x[0])/(2*V[i]);
		B[i][2][5] = (y[0]-y[1])/(2*V[i]);
	}
	static double B_T[2][3*2][3];
	for (int k = 0; k < 2; k++){
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 3; j++) {
				B_T[k][i][j] = B[k][j][i];
			}
		}
	}

	for (int i = 0; i < 4; i++){
		for (int j = 0; j < 4; j++){
			K[i][j] = 0;
		}
	}

	double sum = 0;
	for (int i = 0; i < 2; i++) {
		sum=0;
		double K1[6][3];// K1=B'C= X*k2 * k2*Y
		for (int i2 = 0; i2 < 6; i2++) { //X
			for (int j2 = 0; j2 < 3; j2++) { //Y
				for (int k2 = 0; k2 < 3; k2++) { //X*Y= X*k2 * k2*Y
					sum += B_T[i][i2][k2] * C[k2][j2];
				}
				K1[i2][j2] = sum;
				sum=0;
			}
		}
		double K_ele[6][6];//K_ele= K1*B= B'CB
		sum = 0;
		for (int i2 = 0; i2 < 6; i2++) {
			for (int j2 = 0; j2 < 6; j2++) {
				for (int k2 = 0; k2 < 3; k2++) {
					sum += K1[i2][k2] * B[i][k2][j2] * V[i];
				}
				K_ele[i2][j2] = sum;
				sum = 0;
			}
		}

		//Assembly of local K into the global K
		for (int n = 0; n < 2*3; n++) {
			for (int n1 = 0; n1 < 2*3; n1++) {
				if (ID_ele[i][n] != -1 && ID_ele[i][n1] != -1)
					K[ID_ele[i][n]][ID_ele[i][n1]] += K_ele[n][n1];
			}
		}
	}

	//Sparsize K ********************************************/
	int data_num2 = 0;
	for (int i = 0; i < 4; i++){
		for (int j = 0; j <4; j++){
			if (K[i][j] != 0){
				data_num2 += 1;
			}
		}
	}
	const int data_num=12;//
	static double sparse_K[data_num];
	int row_sparse_K[3];
	static int col_sparse_K[data_num];//
	row_sparse_K[0] = 0;
	counter = 0;
	for (int i = 0; i < 4; i++){
		for (int j = 0; j < 4; j++){
			if (K[i][j] != 0){
				sparse_K[counter] = K[i][j];
				col_sparse_K[counter] = j;
				counter += 1;
			}
			row_sparse_K[i + 1] = counter;
		}
	}
	cout << counter << "\n";
	//**********************************/

	double F[DOF*(node_num - bound_num)];
	for (int i = 0; i < DOF*(node_num - bound_num); i++){
		F[i] = 0;
	}
	F[56] = -1000;//	F[686] = -1000;

	//Solve********************************

	static double p[length];
	static double r[length];
	static double U[length];
	static double t[length];
	double rho = 0;
	double rhos;
	double rho0;
	double alpha;
	int solved = 0;
	for (int i = 0; i < length; i++)
	{
		p[i] = F[i];
		r[i] = F[i];
		U[i] = 0;
	}
	for (int i = 0; i < length; i++)
	{
		rho += r[i] * r[i];
	}
	cout << "rho    " << rho << "\n";
	for (int j = 0; j < length; j++)
	{
		if (solved == 0)
		{
			double sum = 0;
			for (int i = 0; i < length; i++)
			{
			//	for (int k = RK_sp[i]; k < RK_sp[i + 1]; k++)
				{
				//	sum += (p[CK_sp[k]])*(K_sp[k]);
				}
				t[i] = sum;
				sum = 0;
			}
			cout << "t   " << (t[0]) << "    " << t[1] << "    " << t[2] << "    " << t[3] << "\n";
			double PT = 0;
			for (int i = 0; i < length; i++)
			{
				PT += p[i] * t[i];
			}
			cout << "PT   " << PT << "\n";
			alpha = rho / PT;
			cout << "alpha       " << alpha << "\n";
			for (int i = 0; i < length; i++)
			{
				U[i] += alpha*p[i];
				r[i] -= alpha*t[i];
			}
			cout << "U   " << *(U + 56) << "\n" << "\n";
			//cout << "U   " << *(U) << "   " << *(U + 1) << "   " << *(U + 2) << "   " << *(U + 3) << "\n"<<"\n";
			//cout << "r   " << *(r) << "   " << *(r + 1) << "   " << *(r + 2) << "   " << *(r + 3) << "\n";
			//cout << *(U + 1) << "\n";
			rhos = rho;
			rho = 0;
			for (int i = 0; i < length; i++)
			{
				rho += r[i] * r[i];
			}
			if ((rho / rhos) < 0.2)
			{
				solved = 1;
				cout << "HEEEEEEEEEEY" << "\n";
			}
			cout << "rho    " << rho << "\n";
			for (int i = 0; i < length; i++)
			{
				p[i] = r[i] + (rho / rhos)*p[i];
			}
		}
	}



	cout << "Bye Bye" << "\n";
	//system("pause");
	return 0;
}

void Multi2(double **x, int Ix, int Jx, double **y, int Iy, int Jy, double **z){
	//Rerurns x*y & Jx=Iy
	double sum = 0;

	for (int i = 0; i < Ix; i++){
		for (int j = 0; j < Jy; j++){
			for (int k = 0; k < Iy; k++){
				sum += (*(*(x + i) + k))*(*(*(y + k) + j));
			}
			// cout<<sum<<"\n";
			*(*(z + i) + j) = sum;
			sum = 0;
		}
	}
}

void Multi(double **x, int Ix, int Jx, double **y, int Iy, int Jy){
	//Rerurns x*y & Jx=Iy
	double sum = 0;
	double **R;
	R = new double *[Ix];
	for (int i = 0; i < Ix; i++)
		R[i] = new double[Jy];
	for (int i = 0; i < Ix; i++){
		for (int j = 0; j < Jy; j++){
			for (int k = 0; k < Iy; k++){
				sum += (*(*(x + i) + k))*(*(*(y + k) + j));
			}
			// cout<<sum<<"\n";
			*(*(R + i) + j) = sum;
			sum = 0;
		}
	}
}

void SpVec(double *A_sp, int *RA_sp, int *CA_sp, double *b_vec, int length_vec){
	//returns A*b
	static double C[1416];//define it cause it's not a good idea to return local variable adress
	double sum = 0;
	for (int i = 0; i < length_vec; i++){
		for (int j = *(RA_sp + i); j < *(RA_sp + i + 1); j++){
			sum += *(b_vec + *(CA_sp + j))*(*(A_sp + j));
		}
		//cout<<sum<<"\n";
		C[i] = sum;
		sum = 0;
	}
	//return C;
}

double VecVec(double *a, double *b, int n){
	static double L;
	L = 0;
	for (int i = 0; i < n; i++)
		L += *(a + i)*(*(b + i));
	return L;
}

//***************/
void Solve(double *A_sp, int *RA_sp, int *CA_sp, double *B, int length2){
	const int length = 1416;
	double p[length];
	double r[length];
	static double x[length];
	double *t;
	double rho;
	double rhos;
	double rho0;
	double alpha;
	/*
	for (int i = 0; i<length; i++){
	p[i] = B[i];
	r[i] = B[i];
	}
	rho = VecVec(r, r, length);
	cout << "rho    " << rho << "\n";
	for (int j = 0; j<length; j++)
	{
	//t = SpVec(A_sp, RA_sp, CA_sp, p, 3);
	cout << "t   "<< *(t+1) << "\n";
	alpha = rho / VecVec(p, t, length);
	cout<<"alpha       "<<alpha<<"\n";
	for (int i = 0; i<length; i++)
	{
	x[i] += alpha*p[i];
	r[i] -= alpha*(*(t + i));
	}
	cout << *x << "\n";
	rhos = rho;
	rho = VecVec(r, r, length);
	if ((rho / rhos) < 0.02)
	break;
	cout << "rho    " << rho << "\n";
	for (int i = 0; i<length; i++)
	p[i] = r[i] + (rho / rhos)*p[i];
	}
	return x;
	 */
	//return 0;
}
