#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include "project3_library.h"

using namespace std;

ofstream ofile;
ifstream ifile;

int main(int argc, char* argv[]){
	bool test, printmatrix, printvectors, zproj0;
	bool** stateset;

	int n, p, omega, m, count, nchoosem, zpcount;
	int zproj, counter;
	int i, j, k, q, number;
	int* list;

	char* infile;
	char* outfile;

	double energy, g, d, tolerance;
	double** hamiltonian, **vectors, **eigenvalues;
	tolerance = 1e-10;

	if(argc<9){
		cout << "Bad usage. Enter also 'level_count degeneracy particles level_spacing g bool_printmatrix bool_printvectors bool_0_1' on same line." << endl;
		exit(1);
	}
	else{
		p = atoi(argv[1]);
		omega = atoi(argv[2]);
		n = p*omega;
		m = atoi(argv[3]);
		d = atof(argv[4]);
		g = atof(argv[5]);
		printmatrix = atoi(argv[6]);
		printvectors = atoi(argv[7]);
		zproj0 = atoi(argv[8]);
	}

	choose(nchoosem, n, m);
	if(m%2==0){
		if(omega%2==0){
			choose(zpcount, p*omega/2, m/2);
		}
		else{
			choose(zpcount, p*(omega-1)/2, m/2);
		}
	}
	else if(m%2==1){
		zpcount = 0;
	}

	array_alloc(list,zpcount);
	matrix_alloc(stateset, nchoosem, n);
	matrix_alloc(hamiltonian, nchoosem, nchoosem);
	matrix_alloc(vectors, nchoosem, nchoosem);
	matrix_alloc(eigenvalues, nchoosem, 2);

	for(i=0;i<nchoosem;i++){
		vectors[i][i] = 1.0;
	}
	for(i=0;i<nchoosem;i++){
		eigenvalues[i][0] = 0;
		eigenvalues[i][1] = -1;
	}

	construct_stateset(stateset, count, n, m);

/**************************************
Print states to screen	with
	zprojection of state
	even # of particles per level
	zproj==0 && even # of particles
**************************************/
	cout << endl << "States:" << endl;
	for(i=0;i<nchoosem;i++){
		zedangularmomentum(stateset[i], zproj, omega, n);
		detectpairs(stateset[i], test, omega, n);
		cout << setw(10) << i << setw(10);	
		for(j=0;j<n;j++){
			cout << stateset[i][j];
		}
		cout << setw(10) << zproj << setw(10) << test << setw(10) << (zproj==0&&test==1);
		cout << endl;
	}
	cout << endl;

	if(zproj0==true){
		derive_hamiltonian_matrix_zproj0(hamiltonian, stateset, nchoosem, n, omega, d, g);
	}
	else if(zproj0==false){
		derive_hamiltonian_matrix(hamiltonian, stateset, nchoosem, n, omega, d, g);
	}
	if(printmatrix==true){
		cout << endl;
		cout << "Hamiltonian:" << endl;
		for(i=0;i<nchoosem;i++){
			for(j=0;j<nchoosem;j++){
				cout << setw(9) << hamiltonian[i][j];
			}
			cout << endl;
		}
		cout << endl;
	}

	jacobi_simrot_eigen_solver(hamiltonian, vectors, nchoosem, tolerance, count);
	matrix_diag_sort(hamiltonian, eigenvalues, nchoosem);

/************************
Print energies to screen
************************/

	if(printmatrix==1||printvectors==1){
		cout << "original position & energy" << endl;
		for(i=0;i<nchoosem;i++){
			cout << setw(10) << eigenvalues[i][1];
		}
		cout << endl;
		for(i=0;i<nchoosem;i++){
			cout << setw(10) << eigenvalues[i][0];
		}
		cout << endl;
		cout << endl;
	}
	else if(printmatrix==0&&printvectors==0){
		for(i=nchoosem-1;i>-1;i--){
			cout << setw(10) << eigenvalues[i][0] << setw(10) << eigenvalues[i][1];
			cout << endl;
		}
		cout << "energy & original position" << endl;
	}

/************************
Print vectors to screen
************************/
	if(printvectors==true){
		for(i=0;i<nchoosem;i++){
			for(j=0;j<nchoosem;j++){
				k = eigenvalues[j][1];
				cout << setw(10) << vectors[i][k];
			}
			cout << endl;
		}	
		cout << endl;
	}

	array_delete(list);
	matrix_delete(stateset, nchoosem);
	matrix_delete(hamiltonian, nchoosem);
	matrix_delete(vectors, nchoosem);
	matrix_delete(eigenvalues, nchoosem);

	return 0;

}


