#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>



void write_complex_components(int Nc, Complex *ind_com, Fullmol *bases, ofstream &outfile, int it, Parms plist, std::vector<string> &infofilenames) {
	int i, j;
	int size;
	int p1;
	outfile << "iter: " << it << " Ncomplexes:" <<Nc<<endl;
	for (i = 0; i < Nc; i++) {
	  outfile<<i<<' ';
	  for(j=0;j<plist.Nprotypes;j++){
	    
	    if(ind_com[i].NofEach[j]!=0)
	      outfile <<infofilenames[j].substr(0,3)<<' '<<ind_com[i].NofEach[j]<<' ';
	  }
	  outfile<<endl;
	  
	  
	}
	  
}
