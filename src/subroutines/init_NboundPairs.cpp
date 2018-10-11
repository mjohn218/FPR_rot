#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>



void init_NboundPairs(int *nBoundPairs, std::vector<int> &proPairlist, ofstream &outfile,  Parms plist, std::vector<string> &infofilenames, Protein *wholep) {
	int i, j;
	
	int index;
	int p1, p2;
	outfile << "TIME(us)" << '\t';
	for(p1=0;p1<plist.Nprotypes;p1++){
	  for(j=0;j<wholep[p1].npropart;j++){
	    
	    p2=wholep[p1].propart[j];
	    if(p2>=p1){
	      
	      index=p1*plist.Nprotypes+p2;//only store pair of proteins once, so first index (p1) must be <= second index.
	      proPairlist.push_back(index);
	      
	      outfile <<"'"<<infofilenames[p1].substr(0,3)<<","<<infofilenames[p2].substr(0,3)<<"'"<<'\t';
	      cout <<"Pro pair: "<<" p1: "<<p1<<" p2: "<<p2<<" index: "<<index<<' '<<infofilenames[p1].substr(0,3)<<","<<infofilenames[p2].substr(0,3)<<"'"<<'\t';
	    }
	  }
	}
	outfile<<"Nloops"<<endl;
	cout <<"size of proPairlist: "<<proPairlist.size();
	
}
