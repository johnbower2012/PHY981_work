#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include "project2_library.h"

using namespace std;

ofstream ofile;
ifstream ifile;

void operator_annihilate(bool*& state, bool& test, int i){
	if(state[i]==true){
		state[i]=false;
		test = true;
	}
	else if(state[i]==false){
		test = false;
	}
}
void operator_create(bool*& state, bool& test, int i){
	if(state[i]==false){
		state[i]=true;
		test = true;
	}
	else if(state[i]==true){
		test = false;
	}
}
void operator_hamiltonian_0(bool*& state, double& energy, int n, int omega, double d){
	int i, j;
	bool test;
	energy = 0;

	for(i=0;i<n;i++){
		operator_annihilate(state, test, i);
		if(test==true){
			operator_create(state, test, i);
			energy += (i/omega)*d;
		}
	}
}
void operator_hamiltonian_1(bool*& state, double& energy, int r, int s, int n, int omega, double g){
	int level, partner;
	bool test;
	energy = 0;

	operator_annihilate(state, test, r);
	if(test==true){
		level = r/omega;
		partner = (level+1)*omega - 1 - r%omega;
		operator_annihilate(state, test, partner);
		if(test==true){
			operator_create(state, test, s);
			if(test==true){
				level = s/omega;
				partner = (level+1)*omega - 1 - s%omega;
				operator_create(state, test, partner);
				if(test==true){
					energy = -g;
			}
			}
		}
	}
}


int main(int argc, char* argv[]){
	bool test;
	bool** stateset;

	int n, p, omega, m, count, nchoosem, zpcount;
	int zproj, counter;
	int i, j, k, q;
	int* list;

	char* infile;
	char* outfile;

	double energy, g, d, tolerance;
	double** hamiltonian, **vectors, **eigenvalues;
	tolerance = 1e-10;

	if(argc<8){
		cout << "Bad usage. Enter also 'level_count degeneracy particles level_spacing g infilename outfilename' on same line." << endl;
		exit(1);
	}
	else{
		p = atoi(argv[1]);
		omega = atoi(argv[2]);
		n = p*omega;
		m = atoi(argv[3]);
		d = atof(argv[4]);
		g = atof(argv[5]);
		infile = argv[6];
		outfile = argv[7];
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

	cout << endl;
	cout << setw(10) << count << " states" << endl;
	cout << endl;
/**************************************
Print states to screen	with
	zprojection of state
	even # of particles per level
	zproj==0 && even # of particles
**************************************/
	for(i=0;i<nchoosem;i++){
		zedangularmomentum(stateset[i], zproj, omega, n);
		detectpairs(stateset[i], test, omega, n);
		cout << setw(10);		
		for(j=0;j<n;j++){
			cout << stateset[i][j];
		}
		cout << setw(10) << zproj << setw(10) << test << setw(10) << (zproj==0&&test==1);
		cout << endl;
	}
	cout << endl;

/**************************************
Print list of states which meet
	zproj==0 and even # per level
**************************************/
	construct_zproj_pairs_list(stateset, list, 0, nchoosem, omega, n);	
	for(i=0;i<zpcount;i++){
		cout << setw(10) << list[i];
	}
	cout << endl;
	cout << endl;

/**********************
Write hamiltonian file
**********************/
/*
	ofile.open(outfile);
	for(i=0;i<nchoosem;i++){
		for(j=0;j<nchoosem;j++){
			h_pairs(stateset[i], stateset[j], energy, omega, n, g, d);
			ofile << setw(10) << energy;
		}
		ofile << endl;
	}
	ofile.close();
*/
/**********************
Read hamiltonian to be
diagonalized -------
	for now use above
**********************/
/*
	ifile.open(infile);
	for(i=0;i<nchoosem;i++){
		for(j=0;j<nchoosem;j++){
			ifile >> hamiltonian[i][j];
		}
	}
	ifile.close();

	for(i=0;i<nchoosem;i++){
		for(j=0;j<nchoosem;j++){
			cout << setw(10) << hamiltonian[i][j];
		}
		cout << endl;
	}
	cout << endl;

	jacobi_simrot_eigen_solver(hamiltonian, vectors, nchoosem, tolerance, count);
	matrix_diag_sort(hamiltonian, eigenvalues, nchoosem);
*/
/************************
Print energies to screen
************************/
/*
	for(i=0;i<nchoosem;i++){
		cout << setw(10) << eigenvalues[i][1];
	}
	cout << endl;
	for(i=0;i<nchoosem;i++){
		cout << setw(10) << eigenvalues[i][0];
	}
	cout << endl;
	cout << endl;
*/
/************************
Print vectors to screen
************************/
/*
	for(i=0;i<nchoosem;i++){
		for(j=0;j<nchoosem;j++){
			k = eigenvalues[j][1];
			cout << setw(10) << vectors[i][k];
		}
		cout << endl;
	}	
	cout << endl;
*/
	int z;
	double energysum;
	bool*temp;
	array_alloc(temp, n);

	for(z=0;z<nchoosem;z++){
		for(i=0;i<nchoosem;i++){
			energysum = 0;
			for(j=0;j<n;j++){
				for(k=0;k<n;k++){
					for(q=0;q<n;q++){
						temp[q] = stateset[i][q];
					}
					operator_hamiltonian_1(temp, energy, j, k, n, omega, g);
					if(energy==-g){
						overlap(stateset[z], temp, test, n);
						if(test==true){
							energysum += energy;
						}
					}
				}
			}
			energysum /= 4;
			if(i==z){
				operator_hamiltonian_0(stateset[i], energy, n, omega, d);
				energysum += energy;
			}
			cout << setw(5) << energysum;
		}
		cout << endl;
	}
	cout << endl;
	
	array_delete(temp);

	array_delete(list);
	matrix_delete(stateset, nchoosem);
	matrix_delete(hamiltonian, nchoosem);
	matrix_delete(vectors, nchoosem);
	matrix_delete(eigenvalues, nchoosem);

	return 0;

}


