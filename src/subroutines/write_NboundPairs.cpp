#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>



void write_NboundPairs(int *nBoundPairs, std::vector<int> &proPairlist,  ofstream &outfile, int it, Parms plist){
	int i, j;
	
	int index;
	outfile << it*plist.deltat << '\t';
	for (i = 0; i < proPairlist.size(); i++) {
	  index=proPairlist[i];
	  //outfile<<" i in loop: "<<i<<" index: "<<index<<'\t';
	  outfile<<nBoundPairs[index]<<'\t';
	}
	outfile<<plist.nloop<<endl;
	    
}
