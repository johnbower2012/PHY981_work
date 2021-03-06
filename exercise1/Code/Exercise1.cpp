#include<iostream>
#include<cmath>
#include<iomanip>
#include<armadillo>
#include<fstream>

using namespace std;
using namespace arma;

ofstream ofile;
ifstream ifile;

//Initiate main function. Use input file 'infilename' to output two files, 
//'outfilename' for all but 'driplinename' data.
int main(int argc, char* argv[]){
	char* infilename;
	char* outfilename;
	char* driplinename;
	int i, j, k, n, z, ndata, maxn, maxa, maxz;
	double bedata, error, a1, a2, a3, a4, a, snbe, snld, t2,t3;
	
	a1=15.49; a2=17.23; a3=0.679; a4=22.6;
	ndata=2375; maxz=110; maxn=159; maxa=269;
	mat be(maxz+1,maxn+1), liquiddrop(maxz+1,maxn+1);
//Check for proper command line entry
	if(argc<=3){
		cout << "Bad usage. Include: infilename outfilename driplinename." << endl;
		exit(1);
	}
	else if(argc>3){
		infilename=argv[1];
		outfilename=argv[2];
		driplinename=argv[3];
		ifile.open(infilename);
		ofile.open(outfilename);
	}
//Input data from 'infilename', 'bedata.dat' for our purposes here 
//Note that we do not divide by a yet to prevent unnecessary calculations with separation energies which are not given per nucleon.
	for(i=1; i<=ndata; i++){
		ifile >> z; ifile >> a; ifile >> bedata; ifile >> error;
		ifile >> j; ifile >> n;
		be(z,n) = bedata;
	}
//Calculate liquid drop predictions for binding energy
	for(i=1; i<=maxz; i++){
		for(j=1; j<=maxn; j++){
			a=i+j;
			liquiddrop(i,j)=a1*a-a2*pow(a, 2.0/3.0)-a3*pow(i,2.0)/pow(a,1.0/3.0)-a4*pow((j-i),2.0)/a;
		}
	}
//Proton number. Neutron number. Nucleon number. BE input data. BE liquid drop (ld) volume term. BE ld v+surface area. 
//BE ld v+sa+coulomb. BE ld full prediction. abs(Error) between the two. Separation energies from data. SE from ld. abs(Error) between the two.			
	ofile << "# z" << setw(8) << "n" << setw(8) << "a" << setw(15) << "bedata" << setw(15) << "beldv" << setw(15) << "beldv+sa" << setw(15);
	ofile << "beldv+sa+c" << setw(15) << "beldfull" << setw(15) << "error" << setw(15) << "sn_data" << setw(15) << "sn_liquiddrop" << setw(15) << "error" << endl;
	ofile.precision(8);
//Write to file the calculations. Write the separation energies only if both b(i,j) and be(i,j-1) are non-zero and z=8,9,20,28,50, or 82.
//Then add dripline calculations.
	for(i=1; i<=maxz; i++){
		for(j=1; j<=maxn; j++){
			if(be(i,j)!=0.0){
				a= i + j;
				error = abs(be(i,j)-liquiddrop(i,j))/a;
				t2 = a1 - a2*pow(a,-1.0/3.0);
				t3 = t2 - a3*pow(i,2.0)*pow(a,-4.0/3.0);
				ofile << setiosflags(ios::showpoint);
				ofile << setw(3) << i << setw(8);
				ofile << j << setw(8);
				ofile << setprecision(3) << a << setw(15) << setprecision(8);
				ofile << be(i,j)/a << setw(15);
				ofile << a1 << setw(15);
				ofile << t2 << setw(15);
				ofile << t3 << setw(15);
				ofile << liquiddrop(i,j)/a << setw(15);
				ofile << error << setw(15);
				if((i==8||i==9||i==20||i==28||i==50||i==82)&&(be(i,j-1)!=0.0)){
					snbe = be(i,j)-be(i,j-1);
					ofile << snbe << setw(15);
					snld = liquiddrop(i,j)-liquiddrop(i,j-1);
					ofile << snld << setw(15);
					error = abs(snbe-snld);
					ofile << error;
				}
				ofile << endl;
			}
		}
	}

	ifile.close();
	ofile.close();

//Dripline Calculations. New file for simplicity.
	ofile.open(driplinename);
	ofile << "# z" << setw(15) << "n" << endl;
	for(i=1; i<=120; i++){
		j=1;
		a=i+j;
		while(a1-a2*(pow(a, 2.0/3.0)-pow(a-1,2.0/3.0))-a3*pow(i,2.0)*(pow(a,-1.0/3.0)-pow(a-1,-1.0/3.0))-a4*(pow((j-i),2.0)/a-pow((j-i-1),2.0)/(a-1))>0.0){
			j++;
			a=i+j;
		}
		ofile << setw(15) << i << setw(15) << j << endl;	
	}
	ofile.close();

	return 0;
}
