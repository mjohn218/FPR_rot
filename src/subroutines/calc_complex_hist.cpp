#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>



double calc_complex_hist(int Nc, Complex *ind_com, Fullmol *bases, ofstream &outfile, int it, Parms plist, std::vector<string> &infofilenames, int *Ncopy) {
	int i, j;
	int size;
	int p1;
	double index;//might exceed max integer value.
	double mult[plist.Nprotypes];
	int loc, flag;

	std::vector<double> assemblylist;
	std::vector<int> histogram;
	std::vector<int> complexrep;
	for(j=0;j<plist.Nprotypes;j++){
	  mult[j]=1;
	  for(i=0;i<j;i++)
	    mult[j]=mult[j]*(Ncopy[i]+1);
	  //cout <<"mult factor for type: "<<j<<" is: "<<mult[j]<<endl;
	}
	assemblylist.reserve(plist.Nprotypes);
	complexrep.reserve(plist.Nprotypes);
	histogram.reserve(plist.Nprotypes);
	
	//cout <<"Ncomplexes: "<<Nc<<endl;
	/*Create the first complex type.*/
	i=0;
	index=0;
	for(j=0;j<plist.Nprotypes;j++){
	  index+=ind_com[i].NofEach[j]*mult[j];
	}
	
	//new assembly type.
	
	assemblylist.push_back(index);
	complexrep.push_back(i);//example complex with this composition
	histogram.push_back(1);//created a new assembly type with one so far
	//cout <<"New assembly type: "<<index<<" ffrom complex: "<<i<<" size of assemblylist: "<<assemblylist.size()<<endl;
	  
	for (i = 1; i < Nc; i++) {
	  index=0;
	  for(j=0;j<plist.Nprotypes;j++){
	    index+=ind_com[i].NofEach[j]*mult[j];
	  }
	  flag=0;
	  for(int a=0;a<assemblylist.size();a++){
	    if(index==assemblylist[a]){
	      loc=a;
	      flag=1;
	      a=assemblylist.size();//break from loop, found your index already.
	    }
	  }
	  if(flag==0){
	    //new assembly type.
	    
	    assemblylist.push_back(index);
	    complexrep.push_back(i);//example complex with this composition
	    histogram.push_back(1);//created a new assembly type with one so far
	    //   cout <<"New assembly type: "<<index<<" ffrom complex: "<<i<<" size of assemblylist: "<<assemblylist.size()<<endl;
	  }else{
	    histogram[loc]+=1;
	  }
	  
	}
	//	cout <<" Assembly types: "<<assemblylist.size()<<endl;
	//cout <<" histogram of last one: "<<histogram[assemblylist.size()-1]<<endl;
	/*Calculated histograms of the assemblies, now write them out.*/
	double meanComplexSize=0.0;//This will count mean complex size over all complexes >1 protein
	int numComplexTypes=0;
	int totProteins=0;
	outfile << "iter: " << it << endl;
	for(int a=0;a<assemblylist.size();a++){
	  int c1=complexrep[a];
	  outfile<<histogram[a]<<'\t';
	  totProteins=0.0;
	  for(j=0;j<plist.Nprotypes;j++){
	    if(ind_com[c1].NofEach[j]!=0){
	      outfile <<infofilenames[j].substr(0,3)<<": "<<ind_com[c1].NofEach[j]<<". ";
	      totProteins+=ind_com[c1].NofEach[j];
	    }
	  }
	  outfile<<endl;
	  if(totProteins>1){
	    numComplexTypes+=histogram[a];
	    meanComplexSize+=histogram[a]*totProteins;
	    
	  }
	}
	meanComplexSize=meanComplexSize/(1.0*numComplexTypes);
	return meanComplexSize;
}
