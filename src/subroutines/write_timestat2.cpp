#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void write_timestat2(ofstream &outfile, ofstream &molectypesfile, ofstream &timestatfiletext, Protein *wholep, Fullmol *bases, Parms &plist, int *Ncopy, int iter, Complex *ind_com, double deltat, int Nprotypes, int **molectypes, int &currentnumberofmolectypes) {
	int i, j, p, k, uniqmolexists;
	int s, ptype,numuniqmol, uniqcmpid, existflag, existmolposition, compindex;
	int CurrNumComplexTypes=0;
	int HEIGHT = plist.ntotalcomplex;
	int WIDTH = 2+2*Nprotypes; //max number of different molecules in a complex is Nprotypes

	int **adjmat;
	// Allocate memory for complex info
	adjmat = new int*[HEIGHT];
	for (i = 0; i < HEIGHT; i++)
		adjmat[i] = new int[WIDTH];

	for (i=0;i<HEIGHT;i++){
		for(j=0;j<WIDTH;j++){
			adjmat[i][j]=0;
		}
	}

	int **complexstat; //define and allocate unique complex ids/positions
	complexstat = new int*[MAXNUMCOMPLEXTYPES];
	for (i = 0; i < MAXNUMCOMPLEXTYPES; i++)
		complexstat[i] = new int[2];

	for(i=0;i<MAXNUMCOMPLEXTYPES;i++){
		for(j=0;j<2;j++){
			complexstat[i][j]=0;
		}
	}

	for (i=0;i<HEIGHT;i++){

		existflag = 0;
		if(CurrNumComplexTypes==0){
			complexstat[CurrNumComplexTypes][0]=i;
			CurrNumComplexTypes += 1;
		}

		numuniqmol = 0;
		uniqcmpid = 1;
		//		c=bases[i].mycomplex;
		s=ind_com[i].mysize;

		for(j=0;j<s;j++){
			existflag = 0;
			existmolposition = 0;
			p=ind_com[i].plist[j];
			//			cout<<p<<"\t";
			ptype=bases[p].protype;
			for(k=0;k<numuniqmol;k++){
				if(adjmat[i][2*(k+1)] == ptype){
					existflag = 1;
					existmolposition = k;
				}
			}
			if(existflag==1){
				//				adjmat[i][2*(j+1)] = ptype;
				adjmat[i][2*(existmolposition+1)+1] += 1;
			}else{
				numuniqmol += 1;
				adjmat[i][2*numuniqmol] = ptype;
				adjmat[i][2*numuniqmol+1] += 1;
			}

		}
		adjmat[i][1] = numuniqmol;
		for(k=1;k<2*(numuniqmol)+2;k++){ //calculate unique compid
			uniqcmpid = uniqcmpid*(adjmat[i][k]+1);
			//			cout<<adjmat[i][k]<<"\t";
		}
		//		cout<<endl;
		adjmat[i][0] = uniqcmpid;

		existflag=0;
		for(j=0;j<CurrNumComplexTypes;j++){
			if(adjmat[i][0]==adjmat[complexstat[j][0]][0]){
				existflag = 1;
				break;
			}
		}
		if(existflag==0){
			complexstat[CurrNumComplexTypes][0]=i;
			CurrNumComplexTypes += 1;
		}

	}

	outfile << iter << "\t" <<iter * deltat * 1e-6 << "\t";

	for (i=0;i<CurrNumComplexTypes;i++){

		uniqmolexists = 0;
		compindex = complexstat[i][0];
		numuniqmol = adjmat[compindex][1];

		for(j=0;j<currentnumberofmolectypes;j++){
			if(adjmat[compindex][0]==molectypes[j][0]){
				uniqmolexists = 1;
				break;
			}
		}

		if(uniqmolexists==0){//this is for keeping track of unique types of molecules produced

			for(k=0;k<2*adjmat[compindex][1]+2;k++){
				molectypes[j][k] = adjmat[compindex][k];
			}
			for (k = 1; k < numuniqmol+1; k++) {
					molectypesfile <<"P";
					molectypesfile <<adjmat[compindex][2*k];
					molectypesfile << ":";
					molectypesfile <<adjmat[compindex][2*k+1];
			}
			molectypesfile<<endl;
			currentnumberofmolectypes += 1;
		}

		for(j=0;j<HEIGHT;j++){
			if(adjmat[j][0]==adjmat[compindex][0]){
				complexstat[i][1] += 1;
			}
		}

		outfile<<complexstat[i][1];

		for (j = 0; j < plist.Nprotypes; j++) {
			for (k = 1; k < numuniqmol+1; k++) {
				if(adjmat[compindex][2*k]==j){
					timestatfiletext <<"P";
					timestatfiletext <<adjmat[compindex][2*k];
					timestatfiletext << ":";
					timestatfiletext <<adjmat[compindex][2*k+1];
				}
			}
		}

		outfile<<"\t";
		timestatfiletext<<"\t";
	}
	outfile<<endl;
	timestatfiletext<<endl;
	// De-Allocate memory to prevent memory leak
	for (i = 0; i < HEIGHT; i++)
		delete [] adjmat[i];
	delete [] adjmat;

	for (i = 0; i < MAXNUMCOMPLEXTYPES; i++)
		delete [] complexstat[i];
	delete [] complexstat;
}
