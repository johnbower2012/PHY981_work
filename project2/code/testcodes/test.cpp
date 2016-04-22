#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>

using namespace std;

void array_alloc(int*& a, int n);
void array_delete(int*& a);
void matrix_alloc(bool**& a, int rows, int columns);
void matrix_delete(bool**& a, int rows);
void choose(int& choose, int n, int m);
void construct_stateset(bool**& stateset, int& count, int n, int m);
void zedangularmomentum(bool*& a, int& zproj, int omega, int n);
void detectpairs(bool*& a, bool& test, int omega, int n);
void construct_zproj_pairs_list(bool**& stateset, int*& list, int zed_am, int statecount, int omega, int n);


int main(int argc, char* argv[]){
	int n, p, omega, m, nchoosem, count, zpcount, counter;
	int zproj, pairtest;
	int i, j;
	bool test;
	int* list;
	bool** stateset;

	if(argc<4){
		cout << "Bad usage. Enter also 'level_count degeneracy particles on same line." << endl;
		exit(1);
	}
	else{
		p = atoi(argv[1]);
		omega = atoi(argv[2]);
		n = p*omega;
		m = atoi(argv[3]);
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

	construct_stateset(stateset, count, n, m);
	construct_zproj_pairs_list(stateset, list, 0, nchoosem, omega, n);

	for(i=0;i<zpcount;i++){
		cout << setw(10) << list[i];
	}
	cout << endl;

	array_delete(list);
	matrix_delete(stateset, nchoosem);

	return 0;

}

void array_alloc(int*& a, int n){
	int i;
	a = new int[n];
	for(i=0;i<n;i++){
		a[i] = 0;
	}
}
void array_delete(int*& a){
	delete[] a;
}
void matrix_alloc(bool**& a, int rows, int columns){
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
void matrix_delete(bool**& a, int rows){
	for(int i=0;i<rows;i++){
		delete a[i];
	}
	delete[] a;
}
void choose(int& choose, int n, int m){
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
void construct_stateset(bool**& stateset, int& count, int n, int m){
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
void zedangularmomentum(bool*& a, int& zproj, int omega, int n){
	int j;	
	zproj=0;
	for(j=0;j<n;j++){
		if(a[j]==true){
			zproj += (omega - 1 - 2*(j%omega));
		}
	}
}
void detectpairs(bool*& a, bool& test, int omega, int n){
	int j, pairtest;	
	test = true;
	pairtest = 0;
	for(j=0;j<n;j++){
		if(j%omega==0){
			if(pairtest%2==1){
				test = false;				
				break;
			}
			else{
				pairtest = 0;
			}
		}
		if(a[j]==true){
			pairtest++;
		}
	}
	if(pairtest%2==1){
		test = false;
	}
}
void construct_zproj_pairs_list(bool**& stateset, int*& list, int zed_am, int statecount, int omega, int n){
	bool test;	
	int i, counter, zproj;
	counter = 0;
	for(i=0;i<statecount;i++){
		zedangularmomentum(stateset[i], zproj, omega, n);
		detectpairs(stateset[i], test, omega, n);
		if(zproj==zed_am&&test==true){
			list[counter] = i;
			counter++;
		}
	}
}



