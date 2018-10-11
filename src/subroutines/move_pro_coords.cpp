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

void move_pro_coords(int c1, Complex *ind_com, Fullmol *bases, double *chg1) {
	/*Update the position of proteins in complex c1 by the displacement vector chg1
	 Enforce PBC.
	 */

	int j, mp, k;

	for (j = 0; j < ind_com[c1].mysize; j++) {
		mp = ind_com[c1].plist[j];

		for (k = 0; k < bases[mp].ninterface; k++) {
			bases[mp].x[k] += chg1[0];
			bases[mp].y[k] += chg1[1];
			bases[mp].z[k] += chg1[2];
			//      cout<<bases[mp].z[k]<<endl;
		}
		bases[mp].xcom += chg1[0];
		bases[mp].ycom += chg1[1];
		bases[mp].zcom += chg1[2];
		//    cout<<bases[mp].zcom<<endl;

	}

}
