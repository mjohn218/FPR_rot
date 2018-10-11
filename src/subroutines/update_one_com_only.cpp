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

void update_one_com_only(int c1, Complex *ind_com, Fullmol *bases) {
	/*Update Complex COM of a single proteins*/
	double totalmass = 0;

	int s1;
	int mp;
	int i, j;

	ind_com[c1].xcom = 0;
	ind_com[c1].ycom = 0;
	ind_com[c1].zcom = 0;

	s1 = ind_com[c1].mysize;
	/*Size of complex does not change*/
	for (i = 0; i < s1; i++) {
		mp = ind_com[c1].plist[i];
		totalmass += bases[mp].mass;
		
		ind_com[c1].xcom += bases[mp].xcom * bases[mp].mass;
		ind_com[c1].ycom += bases[mp].ycom * bases[mp].mass;
		ind_com[c1].zcom += bases[mp].zcom * bases[mp].mass;

	}
	ind_com[c1].xcom /= totalmass;
	ind_com[c1].ycom /= totalmass;
	ind_com[c1].zcom /= totalmass;

}
