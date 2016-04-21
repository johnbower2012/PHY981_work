#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>

using namespace std;

ofstream ofile;
ifstream ifile;

inline void array_alloc(bool*& a, int n);
inline void array_delete(bool*& a);
inline void array_alloc(int*& a, int n);
inline void array_delete(int*& a);
inline void array_resize(int*& array, int oldsize, int newsize);
inline void matrix_alloc(bool**& a, int rows, int columns);
inline void matrix_delete(bool**& a, int rows);
inline void matrix_alloc(double**& a, int rows, int columns);
inline void matrix_delete(double**& a, int rows);
inline void factorial(int& factorial, int n);
inline void choose(int& choose, int n, int m);
inline void construct_stateset(bool**& stateset, int& count, int n, int m);
inline void print_stateset(bool**& stateset, int statecount, int n);
inline void create(bool*&, int i);
inline void annihilate(bool*& a, int i);
inline void countfilledstates(bool*& a, int& onecount, int n);
inline void overlap(bool*& a, bool*& b, bool& test, int n);
inline void overlap(bool*& a, bool*& b, int& test, int n);	
inline void h0_pairs(bool*& a, bool*& b, double& energy, int omega, int n, double d);
inline void h1_pairs(bool*& a, bool*& b, double& energy, int n, double g);
inline void h_pairs(bool*& a, bool*& b, double& energy, int omega, int n, double g, double d);
inline void zedangularmomentum(bool*& a, int& zedangmom, int omega, int n);
inline void detectpairs(bool*& a, bool& test, int omega, int n);
inline void count_zedproj_pairs(bool**& stateset, int*& list, int zedangmom, int statecount, int omega, int n){
	int i, j, zproj, count;
	bool test;
	count = 0;
	for(i=0;i<statecount;i++){
		zedangularmomentum(stateset[i], zproj, omega, n);
		detectpairs(stateset[i], test, omega, n);
		cout << zproj << endl;
		if(zproj==zedangmom&&test==true){
			list[count] = i;
			count++;
		}
	}
}

void jacobi_simrot_eigen_solver(double**& A, double**& R, int size, double tolerance, int& count);
void matrix_diag_sort(double**& A, double**&B, int size);


int main(int argc, char* argv[]){
	int n, p, omega, m, count, nchoosem, zpcount;
	int onecount, zproj;
	int i, j, k;
	bool test;
	char* infile;
	char* outfile;
	int* list;
	bool** stateset;
	double** hamiltonian, **vectors, **eigenvalues;
	double energy, g, d, tolerance;
	tolerance = 1e-10;
	d = 1;

	if(argc<7){
		cout << "Bad usage. Enter also 'level_count degeneracy particles g infilename outfilename' on same line." << endl;
		exit(1);
	}
	else{
		p = atoi(argv[1]);
		omega = atoi(argv[2]);
		n = p*omega;
		m = atoi(argv[3]);
		g = atof(argv[4]);
		infile = argv[5];
		outfile = argv[6];
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
	cout << setw(10) << nchoosem << " nchoosem" << endl;
	cout << setw(10) << count << " states" << endl;
	cout << endl;

	count_zedproj_pairs(stateset,list,0,nchoosem,omega,n);
	for(i=0;i<nchoosem;i++){
		zedangularmomentum(stateset[i], zproj, omega, n);
		detectpairs(stateset[i], test, omega, n);
		cout << setw(10);
		for(j=0;j<n;j++){
			cout << stateset[i][j];
		}
		cout << setw(10) << zproj << setw(10) << test;
		if(zproj==0&&test==true){ 
			cout << setw(10) << 1;
		}
		else{
			cout << setw(10) << 0;
		}
		cout << endl;
	}
	cout << endl;

	cout << setw(10) << zpcount;
	for(i=0;i<zpcount;i++){
		cout << setw(10) << list[i];
	}
	cout << endl << endl;

	ofile.open(outfile);
	for(i=0;i<nchoosem;i++){
		for(j=0;j<nchoosem;j++){
			h_pairs(stateset[i], stateset[j], energy, omega, n, g, d);
			ofile << setw(10) << energy;
		}
		ofile << endl;
	}
	ofile.close();

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

	for(i=0;i<nchoosem;i++){
		cout << setw(10) << eigenvalues[i][1];
	}
	cout << endl;
	for(i=0;i<nchoosem;i++){
		cout << setw(10) << eigenvalues[i][0];
	}
	cout << endl;
	for(i=0;i<nchoosem;i++){
		for(j=0;j<nchoosem;j++){
			k = eigenvalues[j][1];
			cout << setw(10) << vectors[i][k];
		}
		cout << endl;
	}	
	cout << endl;


	array_delete(list);
	matrix_delete(stateset, nchoosem);
	matrix_delete(hamiltonian, nchoosem);
	matrix_delete(vectors, nchoosem);
	matrix_delete(eigenvalues, nchoosem);

	return 0;

}

inline void array_alloc(bool*& a, int n){
	int i;
	a = new bool[n];
	for(i=0;i<n;i++){
		a[i] = false;
	}
}
inline void array_delete(bool*& a){
	delete[] a;
}
inline void array_alloc(int*& a, int n){
	int i;
	a = new int[n];
	for(i=0;i<n;i++){
		a[i] = 0;
	}
}
inline void array_delete(int*& a){
	delete[] a;
}
inline void array_resize(int*& array, int oldsize, int newsize){
	int* temp = new int[newsize];
	for(int i=0;i<oldsize;i++){
		temp[i]=array[i];
	}

	delete[] array;
	array = temp;
}
inline void matrix_alloc(bool**& a, int rows, int columns){
	int i, j;	
	a = new bool*[rows];
	for(i=0;i<rows;i++){
		a[i] = new bool[columns];
	}
	for(i=0;i<rows;i++){
		for(j=0;j<columns;j++){
			a[i][j] = false;
		}
	}
}
inline void matrix_delete(bool**& a, int rows){
	for(int i=0;i<rows;i++){
		delete a[i];
	}
	delete[] a;
}
inline void matrix_alloc(double**& a, int rows, int columns){
	int i, j;	
	a = new double*[rows];
	for(i=0;i<rows;i++){
		a[i] = new double[columns];
	}
	for(i=0;i<rows;i++){
		for(j=0;j<columns;j++){
			a[i][j] = 0.0;
		}
	}
}
inline void matrix_delete(double**& a, int rows){
	for(int i=0;i<rows;i++){
		delete a[i];
	}
	delete[] a;
}
inline void factorial(int& factorial, int n){
	int i;
	if(n<0){
		cout << "Bad usage of 'factorial.' Must choose a positive number." << endl;
		exit(1);
	}
	factorial = 1;
	for(i=1;i<n+1;i++){
		factorial *= i;
	}
}
inline void choose(int& choose, int n, int m){
	int i, j, k;
	choose = n;
	if(m>n||m<0){
		choose = 0;
	}
	else if(m==n||m==0){
		choose = 1;
	}
	else{
		if(m*2 > n){
			m = n-m;
		}
		for(i=1;i<m;i++){
			choose *= (n-i);
			choose /= (i+1);
		}
	}
}
inline void construct_stateset(bool**& stateset, int& count, int n, int m){
	int i, j;	
	int* m_list;
	
	i=0; count = 0;
	array_alloc(m_list, m);

	for(j=0;j<m;j++){
		m_list[j] = j;
	}

	while(i<m){
		if(m_list[m-1-i]<(n-i)){
			if(i==0){
				for(j=0;j<m;j++){
					stateset[count][m_list[j]] = true;
				}
				m_list[m-1] += 1;
				count++;
			}
			else if(i<m){
				m_list[m-1-i] += 1;
				while(i>0){
					m_list[m-i] = m_list[m-1-i] + 1;
					i -= 1;
				}
			}
		}
		else if(m_list[m-1-i]==n-i){
			i += 1;
		}
	}

	array_delete(m_list);
}
inline void print_stateset(bool**& stateset, int statecount, int n){
	int i, j;
	cout << endl;
	for(i=0;i<statecount;i++){
		for(j=0;j<n;j++){
			cout << stateset[i][j];
		}
		cout << endl;
	}
	cout << endl;
}
inline void create(bool*& a, int i){
	a[i] = a[i]^true;
}
inline void annihilate(bool*& a, int i){
	a[i] = a[i]^true;
}
inline void countfilledstates(bool*& a, int& onecount, int n){	
	onecount = 0;
	int i;
	for(i=0;i<n;i++){
		if(a[i]==true){
			onecount++;
		}
	}
}
inline void overlap(bool*& a, bool*& b, bool& test, int n){
	int stop, i;
	bool* c;
	array_alloc(c,n);

	i=0; stop=0;
	while(i<n&&stop==0){
		c[i] = a[i]^b[i];
		stop += c[i];
		i++;
	}
	if(stop==0){
		test=true;
	}
	else{
		test = false;
	}	

	array_delete(c);
}
inline void overlap(bool*& a, bool*& b, int& test, int n){
	int i;
	bool* c;
	array_alloc(c,n);
	test = 0;

	for(i=0;i<n;i++){
		c[i] = a[i]^b[i];
		test += c[i];
	}

	array_delete(c);
}
inline void h0_pairs(bool*& a, bool*& b, double& energy, int omega, int n, double d){
	bool test;
	int i, level;
	energy = 0;
	overlap(a,b,test,n);
	if(test==true){
		for(i=0;i<n;i++){
			if(a[i]==true){
				level = i/omega;
				energy += 2.0*level*d;
			}
		}
	}
}
inline void h1_pairs(bool*& a, bool*& b, double& energy, int omega, int n, double g){
	int i, test;
	energy = 0;
	overlap(a,b,test,n);
	if(test==0){
		for(i=0;i<n;i++){
			if(a[i]==true){
				energy += -g;
			}
		}
	}
	else if(test==2){
		energy = -g;
	}
}
inline void h_pairs(bool*& a, bool*& b, double& energy, int omega, int n, double g, double d){
	int i, test, level;
	energy = 0;
	overlap(a,b,test,n);
	if(test==0){
		for(i=0;i<n;i++){
			if(a[i]==true){
				level = i/omega;
				energy += 2.0*level*d;
			}
		}
		for(i=0;i<n;i++){
			if(a[i]==true){
				energy += -g;
			}
		}
	}
	if(test==2){
		energy = -g;
	}
}
inline void zedangularmomentum(bool*& a, int& zedangmom, int omega, int n){
	int i, zproj;
	zedangmom = 0;
	for(i=0;i<n;i++){
		if(a[i]==true){
			zproj = i%omega;
			zproj = (omega - 1 - 2*zproj);
			zedangmom += zproj;
		}
	}
}
inline void detectpairs(bool*& a, bool& test, int omega, int n){
	int count, mark, i;
	test = true;
	for(i=0;i<n;i++){
		mark = i%omega;
		if(mark==0){
			if(count%2==1){
				test = false;				
				break;
			}
			count = 0;
		}
		if(a[i]==true){
			count++;
		}
	}
	if(count%2==1){
		test = false;
	}
}
void jacobi_simrot_eigen_solver(double**& A, double**& R, int size, double tolerance, int& count){
	double a_ik, a_il, a_kk, a_ll, r_ik, r_il;
	double tau, tan, cos, sin;
	double temp, A_max;
	A_max = 0.0;
	int i, j, k, l;	
	//create array to find max element of each row in A and its column
	double** a_max = new double*[size];
	for(i=0;i<size;i++){
		a_max[i] = new double[2];
	}
	//initialize array so that any reference to improper value returns error
	for(i=0;i<size;i++){
		a_max[i][0] = 0;
		a_max[i][1] = -1;
	}
	//find max offdiagonal element of each row in A
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
	//set count of similarity rotations to zero for proper count
	count = 0;
	//conduct rotations until maximum offdiagonal in A is below tolerance
	while(A_max>tolerance){
		//calculate cos(theta) and sin(theta) from given offdiagonal
		tau = (A[l][l] - A[k][k])/(2*A[k][l]);
		if(tau>=0) tan = -tau + sqrt(tau*tau+1);
		else if(tau<0) tan = -tau - sqrt(tau*tau+1);
		cos = 1.0/sqrt(tan*tan + 1.0);
		sin = tan*cos;
		//calculate changes to A and R from rotation
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

		//account for change to maximum element and to rows k and l.
		A_max = 0.0;
		a_max[k][0]=0.0;
		a_max[k][1]=-1;
		a_max[l][0]=0.0;
		a_max[l][1]=-1;

		//find new maximum elements for rows k and l
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
		//compare current maximum elements with changes in columns k and l
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
		//find new maximum from array of maximum elements		
		for(i=0;i<size;i++){
			if(a_max[i][0]>A_max){
				A_max = a_max[i][0];
				k = i;
				l = a_max[i][1];
			}
		}
	//+1 similarity rotation
	count++;
	}
	//delete maximum array
	for(i=0;i<size;i++){
		delete[] a_max[i];
	}
	delete[] a_max;
}
//sort diagonal elements of A into given 'size' by '2' matrix B, smallest to largest
void matrix_diag_sort(double**& A, double**&B, int size){
	int i, j, q;
	B[0][0]=A[0][0];
	B[0][1] = 0;
	for(i=1;i<size;i++){
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



