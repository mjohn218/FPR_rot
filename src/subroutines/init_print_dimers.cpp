#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>



void init_print_dimers( ofstream &outfile, int it, Parms plist, std::vector<string> &infofilenames) {
	int i, j;
	int size;
	int p1;
	double index;//might exceed max integer value.
	outfile<<" ITER: "<< "TIME (us) "<<'\t';
	for(j=0;j<plist.Nprotypes;j++){

	  outfile <<"MONO:"<<infofilenames[j].substr(0,3)<<'\t'<<"DIMERS W:"<<infofilenames[j].substr(0,3)<<'\t';
	}
	outfile<<endl;	

}
