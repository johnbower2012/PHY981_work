#include<iostream>
#include<bitset>

using namespace std;

int main(int argc, char* argv[]){
	int n, m, l, i, j, count;

	bitset<16>* stateset;
	count = 0;
	n = atoi(argv[1]);
	m = atoi(argv[2]);
	l = atoi(argv[3]);
	bitset<16> temp;
	cout << temp << endl;
	stateset = new bitset<16> [l];


	for(j=0;j<m;j++){
		temp.set(j,1);
	}
	cout << temp << endl << endl;

	for(i=0;i<(n-m)+1;i++){
		cout << temp << endl;
		stateset[count] = temp;
		count++;
		for(j=m+i;j<n;j++){
			temp.set(j-1,0);
			temp.set(j,1);
			cout << temp << endl;
			stateset[count] = temp;
			count++;
		}
		temp.set(m+i-1,1);
		temp.set(n-1,0);
		temp = temp<<1;
	}
	cout << count << endl << endl;

	for(i=0;i<count;i++){
		cout <<	stateset[i] << endl;
	}

	

	for(i=0;
	delete[] stateset;
	
	return 0;
}
