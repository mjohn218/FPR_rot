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

void set_movestat_zero(int p1, Fullmol *bases, Complex *ind_com, int *movestat) {
	int k1 = bases[p1].mycomplex;
	int i;
	int s1 = ind_com[k1].mysize;
	int mp;

	for (i = 0; i < s1; i++) {
		mp = ind_com[k1].plist[i];
		movestat[mp] = 2; //physically updated position of these proteins
	}

}
