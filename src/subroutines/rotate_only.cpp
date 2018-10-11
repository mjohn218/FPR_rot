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

void rotate_only(int p1, int c1, Complex *ind_com, Fullmol *bases, double *M) {
	/*First rotate all points in complex 1 around the com of protein p1, and translate them by vector d*/
	double *v = new double[3];
	double *v2 = new double[3];

	int mp;
	int i, k;
	int s1 = ind_com[c1].mysize;
	double pivx = bases[p1].xcom;
	double pivy = bases[p1].ycom;
	double pivz = bases[p1].zcom;

	for (i = 0; i < s1; i++) {
		mp = ind_com[c1].plist[i];
		for (k = 0; k < bases[mp].ninterface; k++) {
			v[0] = bases[mp].x[k] - pivx;
			v[1] = bases[mp].y[k] - pivy;
			v[2] = bases[mp].z[k] - pivz;

			rotate(v, M, v2); //includes the interface that will align

			bases[mp].x[k] = pivx + v2[0];
			bases[mp].y[k] = pivy + v2[1];
			bases[mp].z[k] = pivz + v2[2];
		}
		//rotate COM
		v[0] = bases[mp].xcom - pivx;
		v[1] = bases[mp].ycom - pivy;
		v[2] = bases[mp].zcom - pivz;
		rotate(v, M, v2); //includes the interface that will align

		bases[mp].xcom = pivx + v2[0];
		bases[mp].ycom = pivy + v2[1];
		bases[mp].zcom = pivz + v2[2];

	}
	delete[] v;
	delete[] v2;

}
