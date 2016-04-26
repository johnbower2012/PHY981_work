#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include "project3_library.h"

using namespace std;

ofstream ofile;

int main(int argc, char* argv[]){
	bool test, printmatrix, printstates, write;
	bool** stateset, **temp;

	int n, p, omega, m, count, nchoosem, counter;
	int i, j, k, q;
	int* list;

	char* infile;
	char* outfile;

	double energy, g, d, tolerance;
	double** hamiltonian, **hamil, **vectors, **eigenvalues;
	tolerance = 1e-10;

	if(argc<9){
		cout << "Bad usage. Enter also 'level_count degeneracy level_spacing g bool_printstates bool_printmatrix write_print_0_1 outfilename' on same line." << endl;
		exit(1);
	}
	else{
		p = atoi(argv[1]);
		omega = atoi(argv[2]);
		n = p*omega;
		d = atof(argv[3]);
		g = atof(argv[4]);
		printstates = atoi(argv[5]);
		printmatrix = atoi(argv[6]);
		write = atoi(argv[7]);
		outfile = argv[8];
	}
	


	if(write==false){
		cout << endl;
		cout << "p=" << p << setw(10) << "omega=" << omega<< setw(10) << "d=" << d << setw(10) << "g=" << g << endl;
		cout << endl;
		for(i=1;i<n+1;i++){
			m=i;
			choose(nchoosem, n, m);
			matrix_alloc(stateset, nchoosem, n);
			construct_stateset(stateset, count, n, m);
			cout << m << " particles" << endl;
			for(q=0;q<m/2+1;q++){
				count_brokenpairstates(count, n, omega, m, q);
				if(count==0){
					break;
				}

				array_alloc(list,count);
				matrix_alloc(temp, count, n);
				matrix_alloc(hamiltonian, count, count);
				matrix_alloc(hamil, count, count);
				matrix_alloc(vectors, count, count);
				matrix_alloc(eigenvalues, count, 2);
				for(j=0;j<count;j++){
					vectors[j][j] = 1.0;
				}
				for(j=0;j<count;j++){
					eigenvalues[j][0] = 0;
					eigenvalues[j][1] = -1;
				}

				construct_brokenpairs_list(stateset, list, nchoosem, n, omega, q);
				for(j=0;j<count;j++){
					for(k=0;k<n;k++){
						temp[j][k] = stateset[list[j]][k];
					}
				}

				derive_hamiltonian_matrix(hamiltonian, temp, count, n, omega, d, g);
				if(printmatrix==true){
					for(j=0;j<count;j++){
						for(k=0;k<count;k++){
							hamil[j][k] = hamiltonian[j][k];
						}
					}
				}

				jacobi_simrot_eigen_solver(hamiltonian, vectors, count, tolerance, counter);
				matrix_diag_sort(hamiltonian, eigenvalues, count);

				cout << setw(10) << q;
				if(q==1){
					cout << " broken pair";
				}
				else{
					cout << " broken pairs";
				}
				cout << setw(10) << count << " states: " << setw(10) << "energy" << setw(10) << eigenvalues[0][0] << endl;

				if(printstates==true){
					cout << endl;
					for(j=0;j<count;j++){
						cout << setw(10);
						for(k=0;k<n;k++){
							cout << temp[j][k];
						}
						cout << endl;
					}
				}
				if(printmatrix==true){
					cout << endl;
					for(j=0;j<count;j++){
						for(k=0;k<count;k++){
							cout << setw(10) << hamil[j][k];
						}
						cout << endl;
					}
				}
				cout << endl;

				array_delete(list);
				matrix_delete(temp, count);
				matrix_delete(hamil, count);
				matrix_delete(hamiltonian, count);
				matrix_delete(vectors, count);
				matrix_delete(eigenvalues, count);	
			}
			matrix_delete(stateset, nchoosem);
		}
	}
	else if(write==true){
		ofile.open(outfile);
		ofile << "#m";
		for(i=0;i<p/2+1;i++){
			ofile << setw(8) << i << "bp";
		}
		ofile << endl;
		for(i=1;i<n+1;i++){
			m=i;
			choose(nchoosem, n, m);
			matrix_alloc(stateset, nchoosem, n);
			construct_stateset(stateset, count, n, m);
			ofile << m;
			for(q=0;q<m/2+1;q++){
				count_brokenpairstates(count, n, omega, m, q);
				if(count==0){
					break;
				}

				array_alloc(list,count);
				matrix_alloc(temp, count, n);
				matrix_alloc(hamiltonian, count, count);
				matrix_alloc(hamil, count, count);
				matrix_alloc(vectors, count, count);
				matrix_alloc(eigenvalues, count, 2);
				for(j=0;j<count;j++){
					vectors[j][j] = 1.0;
				}
				for(j=0;j<count;j++){
					eigenvalues[j][0] = 0;
					eigenvalues[j][1] = -1;
				}

				construct_brokenpairs_list(stateset, list, nchoosem, n, omega, q);
				for(j=0;j<count;j++){
					for(k=0;k<n;k++){
						temp[j][k] = stateset[list[j]][k];
					}
				}

				derive_hamiltonian_matrix(hamiltonian, temp, count, n, omega, d, g);
				jacobi_simrot_eigen_solver(hamiltonian, vectors, count, tolerance, counter);
				matrix_diag_sort(hamiltonian, eigenvalues, count);

				ofile << setw(10) << eigenvalues[0][0];

				array_delete(list);
				matrix_delete(temp, count);
				matrix_delete(hamiltonian, count);
				matrix_delete(vectors, count);
				matrix_delete(eigenvalues, count);	
			}

			ofile << endl;
		}
	
		matrix_delete(stateset, nchoosem);
	}

	ofile.close();

	return 0;

}


