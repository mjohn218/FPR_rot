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

void put_back_inside(int p1, Parms plist, Fullmol *bases, int *npb, int Ncell, Complex *ind_com, int iter, int &loopit) {

  int j;
  int mybin;
	int xbin, ybin, zbin;
	double small = 1E-6;
	int k1, s1, f;
	double chgx, chgy, chgz;
	double pdx, pdy, pdz;

	int mp;
		
	k1 = bases[p1].mycomplex;
	cout << "Complex: " << k1 << ' ' << ind_com[k1].xcom << ' ' << ind_com[k1].ycom << ' ' << ind_com[k1].zcom << endl;
	cout << "Complex size and radius: " << ind_com[k1].mysize << ' ' << ind_com[k1].radR << endl;
	//write_crds_complex(k1, ind_com, bases);
	//So put back inside
	chgx = 0;
	chgy = 0;
	chgz = 0;
	mybin = 0;
	pdx = bases[p1].xcom - plist.xboxl / 2.0;
	if (pdx > 0) {
	  chgx = -pdx - small;
	} else if (pdx < -plist.xboxl) {
	  chgx = -(pdx + plist.xboxl) + small;
	}
	pdy = bases[p1].ycom - plist.yboxl / 2.0;
	if (pdy > 0) {
	  chgy = -pdy - small;
	} else if (pdy < -plist.yboxl) {
	  chgy = -(pdy + plist.yboxl) + small;
	}
	pdz = bases[p1].zcom - plist.zboxl / 2.0;
	if (pdz > 0) {
	  chgz = -pdz - small;
	} else if (pdz < -plist.zboxl) {
	  chgz = -(pdz + plist.zboxl) + small;
	}
	s1 = ind_com[k1].mysize;
	ind_com[k1].xcom += chgx;
	ind_com[k1].ycom += chgy;
	ind_com[k1].zcom += chgz;
	
	for (j = 0; j < s1; j++) {
	  mp = ind_com[k1].plist[j];
	  
	  bases[mp].xcom += chgx;
	  bases[mp].ycom += chgy;
	  bases[mp].zcom += chgz;
	  
	  //update interface coords
	  for (f = 0; f < bases[mp].ninterface; f++) {
	    bases[mp].x[f] += chgx;
	    bases[mp].y[f] += chgy;
	    bases[mp].z[f] += chgz;
	  }
	}
	cout << "changed positions by: : " << chgx << ' ' << chgy << ' ' << chgz << endl;

	
	loopit++;
	if (loopit > 1000) {
	  cout << "can't fit in box !" << endl;
	  exit(1);
	}

}
