#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>



void print_dimers(int Nc, Complex *ind_com, Fullmol *bases, ofstream &outfile, int it, Parms plist, std::vector<string> &infofilenames, int *Ncopy) {
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
	
	int monomers[plist.Nprotypes];
	int dimers[plist.Nprotypes];
	for(int a=0;a<plist.Nprotypes;a++){
	  monomers[a]=0;
	  dimers[a]=0;
	}
	for(int a=0;a<assemblylist.size();a++){
	  int c1=complexrep[a];
	  if(ind_com[c1].mysize==1){
	    /*This is a monomer*/
	    //outfile<<histogram[a]<<'\t';
	    for(j=0;j<plist.Nprotypes;j++){
	      if(ind_com[c1].NofEach[j]==1){
		monomers[j]=histogram[a];
	      }
	    }
	  }else if(ind_com[c1].mysize==2){
	    /*This is a dimer, either homo (2 copies) or hetero (1 copy)*/
	    for(j=0;j<plist.Nprotypes;j++){
	      if(ind_com[c1].NofEach[j]==1 || ind_com[c1].NofEach[j]==2){
		dimers[j]+=histogram[a];//SUM OVER ALL POSSIBLE DIMERS FOR PROTEIN J
	      }
	    }
	  }
	      //outfile <<infofilenames[j].substr(0,3)<<": "<<ind_com[c1].NofEach[j]<<". ";
	  
	  //	  outfile<<endl;
	    
	  
	}//loop over all assemblies in the system
	outfile  << it << '\t'<<it*plist.deltat<<'\t';  
	for(j=0;j<plist.Nprotypes;j++){
	  outfile<<monomers[j]<<'\t'<<dimers[j]<<'\t';//endl;
	}
	outfile<<endl;
}
