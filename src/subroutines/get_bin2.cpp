#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include <sys/time.h>
#include "reactions.h"

#include "rand_gsl.h"
#include "Faddeeva.hh"
#include "utility_calls.h"
#include "vector_rot_calls.h"

void get_bin2(Parms plist, Fullmol *bases, double cellx, double celly, double cellz, int Nx, int Ny, int Nz, std::vector< std::vector<int> > &binlist, int *npb, int Ncell, Complex *ind_com, int iter) {
	int i, j;
	int mybin;
	int xbin, ybin, zbin;
	double small = 1E-6;
	int k1, s1, f;
	double chgx, chgy, chgz;
	double pdx, pdy, pdz;
	int loopit = 0;
	int mp;
	for (j = 0; j < Ncell; j++){
	  npb[j] = 0;
	  binlist[j].clear();
	}
	
	for (i = 0; i < plist.Ntotalmol; i++) {
	  xbin = int((bases[i].xcom + plist.xboxl / 2.0) / cellx);
	  ybin = int((bases[i].ycom + plist.yboxl / 2.0) / celly);
	  zbin = int(-(bases[i].zcom + small - plist.zboxl / 2.0) / cellz);
	  
	  mybin = xbin + ybin * Nx + zbin * Nx * Ny;
	  //cout <<" protein: "<<i<<" Dz:" <<bases[i].Dz<<" Bins: "<<mybin<<" x, y, z: "<<xbin<<' '<<ybin<<' '<<zbin<<endl;
	  if(bases[i].Dz==0){
	    if(bases[i].zcom-0.1>-plist.zboxl/2.0){
	      cout <<"Protein: "<<i<<" with Dz=0 is off the membrane: "<<bases[i].zcom<<endl;
	      write_crds_complex(bases[i].mycomplex, ind_com, bases);
	      exit(1);
	    }
	  }
	  if(bases[i].zcom>plist.zboxl/2.0 || bases[i].zcom+small<-plist.zboxl/2.0){
	    cout << "OUTSIDE ARRAY, Z: " << i << " " << bases[i].xcom << ' ' << bases[i].ycom << ' ' << bases[i].zcom <<' '<< iter<<endl;
	    put_back_inside(i, plist, bases, npb, Ncell, ind_com, iter, loopit);
	    i=0;//restart bin assignments.
	    for (j = 0; j < Ncell; j++){
	      npb[j] = 0;
	      binlist[j].clear();
	    }
	  }else if(bases[i].ycom>plist.yboxl/2.0 || bases[i].ycom+small<-plist.yboxl/2.0){
	    cout << "OUTSIDE ARRAY, Y: " << i << " " << bases[i].xcom << ' ' << bases[i].ycom << ' ' << bases[i].zcom <<' '<< iter<<endl;
	    put_back_inside(i, plist, bases, npb, Ncell, ind_com, iter, loopit);
	    i=0;//restart bin assignments.
	    for (j = 0; j < Ncell; j++){
	      npb[j] = 0;
	      binlist[j].clear();
	    }
	  }else if(bases[i].xcom>plist.xboxl/2.0 || bases[i].xcom+small<-plist.xboxl/2.0){
	    cout << "OUTSIDE ARRAY, X: " << i << " " << bases[i].xcom << ' ' << bases[i].ycom << ' ' << bases[i].zcom <<' '<< iter<<endl;
	    put_back_inside(i, plist, bases, npb, Ncell, ind_com, iter, loopit);
	    i=0;//restart bin assignments.
	    for (j = 0; j < Ncell; j++){
	      npb[j] = 0;
	      binlist[j].clear();
	    }

	  }else if (mybin > (Ncell - 1) || mybin < 0) {
	    
	    cout << "OUTSIDE ARRAY, mybin: " << i << " " << bases[i].xcom << ' ' << bases[i].ycom << ' ' << bases[i].zcom <<' '<< iter<<endl;
	    put_back_inside(i, plist, bases, npb, Ncell, ind_com, iter, loopit);
	    i=0;//restart bin assignments.
	    for (j = 0; j < Ncell; j++){
	      npb[j] = 0;
	      binlist[j].clear();
	    }
	
	  }else{
	    /*everything is fine.*/
	    bases[i].mybin = mybin;
	    bases[i].mybinind = npb[mybin];
	    binlist[mybin].push_back(i);// * MAXPERBIN + npb[mybin]] = i;
	    npb[mybin]++;//we can replace this with binlist[mybin].size();
	  }
	}

}
