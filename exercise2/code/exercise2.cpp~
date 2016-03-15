#include<iostream>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<armadillo>

using namespace std;
using namespace arma;

int main(){

	double d,g;
	int i,j;
	d=1.0;
	g=-10;
	mat a(3,3); mat r(3,3);
	vec eig(3);
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			if(i==j){
				a(i,j)=2*(i+1)*d-g;
				r(i,j) = 1.0;
			}
			else{
				a(i,j) = -g;
				r(i,j) = 0.0;
			}		
		}
	}
			
	
	
	eig_sym(eig, r, a);

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			cout << a(i,j) << '\t';
		}
		cout << endl;
	}
	cout << endl;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			cout << r(i,j) << '\t';
		}
		cout << endl;
	}
	cout << endl;
	for(i=0;i<3;i++){
		cout << eig(i) << '\t';
	}
	cout << endl;

/*	double** row;
	double* column;
	int n;
	char* infilename;
	char* outfilename;

	if(argc<4){
		cout << "Bad Usage. Input matrix_size infilename outfilename as well on command line." << endl;
		exit(1);
	}
	else{
		n = atoi(argv[1]);
		infilename = argv[2];
		outfilename = argv[3];
	}

	for(int i=1
	row = new (nothrow) double column[n];

	delete[] row;
	delete[] column;
*/

	return 0;
}
