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

void rotate_and_translate_int(int p1, int c1, Complex *ind_com, Fullmol *bases, double *M, double *dtrans, int iind) {

	/*First rotate all points in complex 1 around the 
	 reacting interface of protein p1, and translate them by vector d*/
	double *v = new double[3];
	double *v2 = new double[3];
	double *newpos = new double[3];
	int mp;
	int i, k;
	int s1 = ind_com[c1].mysize;
	double pivx = bases[p1].x[iind];
	double pivy = bases[p1].y[iind];
	double pivz = bases[p1].z[iind];

	for (i = 0; i < s1; i++) {
		mp = ind_com[c1].plist[i];
		for (k = 0; k < bases[mp].ninterface; k++) {
			v[0] = bases[mp].x[k] - pivx;
			v[1] = bases[mp].y[k] - pivy;
			v[2] = bases[mp].z[k] - pivz;

			rotate(v, M, v2); //includes the interface that will align
			translate(v2, dtrans, newpos);
			bases[mp].x[k] = pivx + newpos[0];
			bases[mp].y[k] = pivy + newpos[1];
			bases[mp].z[k] = pivz + newpos[2];
		}
		//rotate COM
		v[0] = bases[mp].xcom - pivx;
		v[1] = bases[mp].ycom - pivy;
		v[2] = bases[mp].zcom - pivz;
		rotate(v, M, v2); //includes the interface that will align
		translate(v2, dtrans, newpos);
		bases[mp].xcom = pivx + newpos[0];
		bases[mp].ycom = pivy + newpos[1];
		bases[mp].zcom = pivz + newpos[2];

	}
	delete[] v;
	delete[] v2;
	delete[] newpos;
}
