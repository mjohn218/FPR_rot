#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "reactions.h"

#include "rand_gsl.h"
#include "Faddeeva.hh"
#include "utility_calls.h"
#include "vector_rot_calls.h"

void move_zcrds(Fullmol *bases, int p1, int p2, double delz1, double delz2, Complex *ind_com) {
	/*move the z coordinates only of complexes for p1 and p2*/
	int i, k;
	int mp;
	int c1 = bases[p1].mycomplex;
	int s1 = ind_com[c1].mysize;
	int c2 = bases[p2].mycomplex;
	int s2 = ind_com[c2].mysize;

	/*translate in z all the components of complex 1 and complex 2*/
	for (i = 0; i < s1; i++) {
	  mp = ind_com[c1].plist[i];
	  for (k = 0; k < bases[mp].ninterface; k++) {
	    bases[mp].z[k] += delz1;
	  }
	  //COM
	  bases[mp].zcom += delz1;
	}
	for (i = 0; i < s2; i++) {
	  mp = ind_com[c2].plist[i];
	  for (k = 0; k < bases[mp].ninterface; k++) {
	    bases[mp].z[k] += delz2;
	  }
	  //COM
	  bases[mp].zcom += delz2;
	}
	
}
