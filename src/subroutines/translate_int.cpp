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

void translate_int(int p1, int c1, Complex *ind_com, Fullmol *bases, double *dtrans) {
	/*First rotate all points in complex 1 around the com of protein p1, and translate them by vector d*/
	double *v = new double[3];
	double *v2 = new double[3];
	double *newpos = new double[3];
	int mp;
	int i, k;
	int s1 = ind_com[c1].mysize;

	for (i = 0; i < s1; i++) {
		mp = ind_com[c1].plist[i];
		for (k = 0; k < bases[mp].ninterface; k++) {
			bases[mp].x[k] += dtrans[0];
			bases[mp].y[k] += dtrans[1];
			bases[mp].z[k] += dtrans[2];
		}

		bases[mp].xcom += dtrans[0];
		bases[mp].ycom += dtrans[1];
		bases[mp].zcom += dtrans[2];

	}
	delete[] v;
	delete[] v2;
	delete[] newpos;
}
