#include<iostream>
#include<iomanip>
#include<cmath>
#include<fstream>
#include "time.h"

using namespace std;

void array_alloc(double*& a, int length);
void matrix_alloc(double**& A, int rows, int columns);
void array_delete(double*& a);
void matrix_delete(double**& A, int rows);
void symmat_offdiag_max(double**& A, int size, double& a_max, int& row, int& column);
void jacobi_simrot_eigen_solver(double**& A, double**& R, int size, double tolerance, int& count);
void array_print(double*& a, int length);
void matrix_print(double**& A, int rows, int columns);
void matrix_print_mult(double**& A, double**& B, int rowsA, int columnsA, int columnsB);
void matrix_diag_sort(double**& A, double**&B, int size);
void init_1eHO_matrices(double**& D, double**& R, int size, double rho_max, double rho_min);
void init_2eHO_matrices(double**& D, double**& R, int size, double rho_max, double rho_min, double omega);
void init_2eHOinter_matrices(double**& D, double**& R, int size, double rho_max, double rho_min, double omega);

ofstream ofile;

int main(int argc, char* argv[]){
	double** D, **R, **eigenvalues;
	char* outfilename;
	double epsilon, V, W, tolerance, time;
	int eigvec_yesno = 1;
	int i, j, k, rows, columns, n, count;
	clock_t start, finish;
	tolerance = 1e-10;
	rows = 5;
	columns = 5;
//Check proper entry arguments
	if(argc<5){
		cout << "Bad usage. Enter also 'outfilename epsilon V W' on same line." << endl;
		exit(1);
	}
	else{
		outfilename = argv[1];
		epsilon = atof(argv[2]);
		V = atof(argv[3]);
		W = atoi(argv[4]);
	}
//allocate matrices
	matrix_alloc(eigenvalues,rows,2);
	matrix_alloc(D, rows, columns);
	matrix_alloc(R, rows, columns);
//Initialize matrices
	for(i=0;i<rows;i++){
		if(i==0||i==4){
			D[i][i] = (i-2)*epsilon;
		}
		else if(i==1||i==3){
			D[i][i] = (i-2)*epsilon + 3.0*W;
		}
		else{
			D[i][i] = (i-2)*epsilon + 4.0*W;
		}
		for(j=0;j<rows;j++){
			if(abs(i-j)==2){
				if(i==2||j==2){
					D[i][j] = sqrt(6.0)*V;
				}
				else if(abs(i)==1||abs(j)==1){
					D[i][j] = 3*V;
				}
			}
			else if(i!=j){
				D[i][j] = 0.0;
			}
			R[i][j] = 0.0;
		}
		R[i][i] = 1.0;
	}
	for(i=0;i<rows;i++){
		eigenvalues[i][0] = 0.0;
		eigenvalues[i][1] = -1.0;
	}
	matrix_print(D,rows,rows);

//Conduct similarity rotations until all offdiag elements < tolerance

	jacobi_simrot_eigen_solver(D, R, rows, tolerance, count);

//Sort the eigenvalues from D from smallest to largest and note original positions
	matrix_diag_sort(D, eigenvalues, rows);
//Write output file
	ofile.open(outfilename);
	ofile.precision(6);
	if(eigvec_yesno==1){
		ofile << "#_i" << setw(15) << "eigenvalues";
		for(i=0;i<rows;i++){
			ofile << setw(14) << "v_" << i;
		}
		ofile << endl;
		for(i=0;i<rows;i++){
			ofile << fixed << i << setw(15) << eigenvalues[i][0] << setw(15) << scientific;
			for(j=0;j<columns;j++){
				k=eigenvalues[j][1];
				ofile << setw(15) << R[i][k];
			}
			ofile << endl;
		}
	}
	else if(eigvec_yesno==0){
		ofile << "#_i" << setw(15) << "eigenvalues" << setw(15) << "j" << endl;
		for(i=0;i<rows;i++){
			ofile << i << setw(15) << eigenvalues[i][0] << setw(15) << eigenvalues[i][1] << endl;
		}
	}
	ofile.close();
	cout.precision(5);
	cout << endl;
	matrix_print(eigenvalues,rows,2);

	matrix_delete(eigenvalues,rows);
	matrix_delete(D,rows);
	matrix_delete(R,rows);

	return 0;
}


/******************************
Begin function definitions
*******************************/


void array_alloc(double*& a, int length){
	a = new double[length];
}
void matrix_alloc(double**& A, int rows, int columns){
	A = new double*[rows];
	for(int i=0;i<rows;i++){
		A[i] = new double[columns];
	}
}
void array_delete(double*& a){
	delete[] a;
}
void matrix_delete(double**& A, int rows){
	for(int i=0;i<rows;i++){
		delete[] A[i];
	}
	delete[] A;
}
void symmat_offdiag_max(double**& A, int size, double& a_max, int& row, int& column){
	double temp = 0.0;
	for(int i=0;i<size;i++){
		for(int j=i+1;j<size;j++){
			temp = fabs(A[i][j]);
			if(temp>a_max){
				a_max = temp;
				row = i;
				column = j;
			}
		}
	}
}
void jacobi_simrot_eigen_solver(double**& A, double**& R, int size, double tolerance, int& count){
	double a_ik, a_il, a_kk, a_ll, r_ik, r_il;
	double tau, tan, cos, sin;
	double temp, A_max;
	A_max = 0.0;
	int i, j, k, l;	
	double** a_max = new double*[size];
	for(i=0;i<size;i++){
		a_max[i] = new double[2];
	}
	for(i=0;i<size;i++){
		a_max[i][0] = 0;
		a_max[i][1] = -1;
	}
	for(i=0;i<size;i++){
		for(j=i+1;j<size;j++){
			temp = fabs(A[i][j]);
			if(temp>a_max[i][0]){
				a_max[i][0] = temp;
				a_max[i][1] = j;
			}
		}
		if(a_max[i][0]>A_max){
			A_max = a_max[i][0];
			k = i;
			l = a_max[i][1];
		}
	}
	count = 0;
	while(A_max>tolerance){
		tau = (A[l][l] - A[k][k])/(2*A[k][l]);
		if(tau>=0) tan = -tau + sqrt(tau*tau+1);
		else if(tau<0) tan = -tau - sqrt(tau*tau+1);
		cos = 1.0/sqrt(tan*tan + 1.0);
		sin = tan*cos;

		a_kk = A[k][k];
		a_ll = A[l][l];
		A[k][k] = a_kk*cos*cos + a_ll*sin*sin - 2.0*sin*cos*A[k][l];
		A[l][l] = a_kk*sin*sin + a_ll*cos*cos + 2.0*sin*cos*A[k][l];
		A[k][l] = 0.0;
		A[l][k] = 0.0;
		for(i=0;i<size;i++){
			if(i!=k && i!=l){
				a_ik = A[i][k];
				a_il = A[i][l];
				A[i][k] = a_ik*cos - a_il*sin;
				A[k][i] = A[i][k];
				A[i][l] = a_il*cos + a_ik*sin;
				A[l][i] = A[i][l];
			}
			r_ik = R[i][k];
			r_il = R[i][l];
			R[i][k] = r_ik*cos - r_il*sin;
			R[i][l] = r_ik*sin + r_il*cos;
		}

	
		A_max = 0.0;
		a_max[k][0]=0.0;
		a_max[k][1]=-1;
		a_max[l][0]=0.0;
		a_max[l][1]=-1;

		for(j=k+1;j<size;j++){
			temp = fabs(A[k][j]);
			if(temp>a_max[k][0]){
				a_max[k][0] = temp;
				a_max[k][1] = j;
			}
		}
		for(j=l+1;j<size;j++){
			temp = fabs(A[l][j]);
			if(temp>a_max[l][0]){
				a_max[l][0] = temp;
				a_max[l][1] = j;
			}
		}
		for(i=0;i<k;i++){
			temp = fabs(A[i][k]);
			if(temp>a_max[i][0]){
				a_max[i][0] = temp;
				a_max[i][1] = k;
			}
		}
		for(i=0;i<l;i++){
			temp = fabs(A[i][l]);
			if(temp>a_max[i][0]){
				a_max[i][0] = temp;
				a_max[i][1] = l;
			}
		}			
		for(i=0;i<size;i++){
			if(a_max[i][0]>A_max){
				A_max = a_max[i][0];
				k = i;
				l = a_max[i][1];
			}
		}
	count++;
	}

	for(i=0;i<size;i++){
		delete[] a_max[i];
	}
	delete[] a_max;
}
void array_print(double*& a, int length){
	for(int i=0;i<length;i++){
		cout << a[i] << endl;
	}
}
void matrix_print(double**& A, int rows, int columns){
	for(int i=0;i<rows;i++){
		for(int j=0;j<columns;j++){
			cout << setw(15);
			cout << A[i][j];
		}
		cout << endl;
	}
}
void matrix_print_mult(double**& A, double**& B, int rowsA, int columnsA, int columnsB){
	double temp;	
	for(int i=0;i<rowsA;i++){
		for(int j=0;j<columnsB;j++){
			temp = 0.0;
			for(int p=0;p<columnsA;p++){
				temp += A[i][p]*B[p][j];
			}
			cout << temp << '\t';
		}
		cout << '\n';
	}
	cout << endl;
}
void matrix_diag_sort(double**& A, double**&B, int size){
	int i, j, q;	
	for(i=0;i<size;i++){
		if(i==0){
			B[0][0]=A[0][0];
			B[0][1] = 0;
		}
		else{
			if(A[i][i]>=B[i-1][0]){
				B[i][0]=A[i][i];
				B[i][1]=i;
			}
			else{
				for(j=0;j<i;j++){
					if(A[i][i]<B[j][0]){
						for(q=i;q>j;q--){
							B[q][0] = B[q-1][0];
							B[q][1] = B[q-1][1];
						}
						B[j][0] = A[i][i];
						B[j][1] = i;
						break;
					}
				}
			}
		}
	}	
}
void init_1eHO_matrices(double**& D, double**& R, int size, double rho_max, double rho_min){
	int i, j;	
	double rho, rho_sq, h, h_sq, temp1, temp2;
	h = (rho_max - rho_min)/((double) (size+1));
	h_sq = h*h;
	temp1 = 2.0/h_sq;
	temp2 = -1.0/h_sq;
	for(i=0;i<size;i++){
		rho = rho_min + (i+1)*h;
		rho_sq = rho*rho;
		for(j=0;j<size;j++){
			if(i==j){
				D[i][j] = temp1 + rho_sq;
				R[i][j] = 1.0;
			}
			else if(i==j+1||i==j-1){
				D[i][j] = temp2;
			}
			else{
				D[i][j] = 0.0;
				R[i][j] = 0.0;
			}
		}
	}	
}
void init_2eHO_matrices(double**& D, double**& R, int size, double rho_max, double rho_min, double omega){
	int i, j;	
	double rho, h, h_sq, temp1, temp2, omega_sq, potential;
	h = (rho_max - rho_min)/((double) (size+1));
	h_sq = h*h;
	omega_sq = omega*omega;
	temp1 = 2.0/h_sq;
	temp2 = -1.0/h_sq;
	for(i=0;i<size;i++){
		rho = rho_min + (i+1)*h;
		potential = omega_sq*rho*rho;
		for(j=0;j<size;j++){
			if(i==j){
				D[i][j] = temp1 + potential;
				R[i][j] = 1.0;
			}
			else if(i==j+1||i==j-1){
				D[i][j] = temp2;
			}
			else{
				D[i][j] = 0.0;
				R[i][j] = 0.0;
			}
		}
	}
}
void init_2eHOinter_matrices(double**& D, double**& R, int size, double rho_max, double rho_min, double omega){
	int i, j;	
	double rho, h, h_sq, temp1, temp2, omega_sq, potential;
	h = (rho_max - rho_min)/((double) (size+1));
	h_sq = h*h;
	omega_sq = omega*omega;
	temp1 = 2.0/h_sq;
	temp2 = -1.0/h_sq;
	for(i=0;i<size;i++){
		rho = rho_min + (i+1)*h;
		potential = omega_sq*rho*rho + 1.0/rho;
		for(j=0;j<size;j++){
			if(i==j){
				D[i][j] = temp1 + potential;
				R[i][j] = 1.0;
			}
			else if(i==j+1||i==j-1){
				D[i][j] = temp2;
			}
			else{
				D[i][j] = 0.0;
				R[i][j] = 0.0;
			}
		}
	}
}
