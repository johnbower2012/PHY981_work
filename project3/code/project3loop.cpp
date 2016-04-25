#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include "project3_library.h"
#include "time.h"

using namespace std;

ofstream ofile;

int main(int argc, char* argv[]){
	bool test, printmatrix, printvectors;
	bool** stateset;

	int n, p, omega, m, count, nchoosem, zpcount;
	int zproj, counter;
	int i, j, k, q;
	int* list;

	char* outfile;

	double energy, g, d, tolerance, time;
	double** hamiltonian, **vectors, **eigenvalues;
	tolerance = 1e-10;

	clock_t start, finish;

	if(argc<7){
		cout << "Bad usage. Enter also 'level_count degeneracy particles level_spacing g outfilename' on same line." << endl;
		exit(1);
	}
	else{
		p = atoi(argv[1]);
		omega = atoi(argv[2]);
		n = p*omega;
		m = atoi(argv[3]);
		d = atof(argv[4]);
		g = atof(argv[5]);
		outfile = argv[6];
	}
	
	ofile.open(outfile);
	ofile << "#p=" << p << setw(10) << "omega=" << omega << setw(10) << "d=" << d << setw(10) << "g=" << g << endl;
	ofile << "#part" << setw(10) << "gs en." << endl;
	cout << "#p=" << p << setw(10) << "omega=" << omega << setw(10) << "d=" << d << setw(10) << "g=" << g << endl;

	for(i=1;i<n+1;i++){
		m = i;
		choose(nchoosem, n, m);
		matrix_alloc(stateset, nchoosem, n);
		matrix_alloc(hamiltonian, nchoosem, nchoosem);
		matrix_alloc(vectors, nchoosem, nchoosem);
		matrix_alloc(eigenvalues, nchoosem, 2);

		start = clock();

		construct_stateset(stateset, count, n, m);
		derive_hamiltonian_matrix(hamiltonian, stateset, nchoosem, n, omega, d, g);

		finish = clock();		

		jacobi_simrot_eigen_solver(hamiltonian, vectors, nchoosem, tolerance, count);
		matrix_diag_sort(hamiltonian, eigenvalues, nchoosem);

		ofile << i << setw(10) << eigenvalues[0][0] << endl;

		time = (finish - start)/((double) CLOCKS_PER_SEC);
		cout << i << " particles" << setw(10) << time << " seconds" << endl;

		matrix_delete(stateset, nchoosem);
		matrix_delete(hamiltonian, nchoosem);
		matrix_delete(vectors, nchoosem);
		matrix_delete(eigenvalues, nchoosem);
	}

	ofile.close();

	return 0;

}


