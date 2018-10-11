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

int measure_overlap(int c1, int c2, Complex *ind_com, Fullmol *bases, int pclath) {

	int t1 = 0;
	int t2 = 0;
	int k;
	int j, mp, mp2;
	int flag = 0;
	int s1 = ind_com[c1].mysize;
	int s2 = ind_com[c2].mysize;
	double tol2 = 1.0;
	double dx, dy, dz;
	double r2;
	if (s1 > 1 && s2 > 1) {

		//    cout <<"Size of complex: "<<s1<<endl;
		//otherwise ther's only 2 proteins in this complex
		for (j = 0; j < s1; j++) {
			mp = ind_com[c1].plist[j];
			if (bases[mp].protype == pclath) {
				//center for clathrin is COM.
				for (k = 0; k < s2; k++) {
					mp2 = ind_com[c2].plist[k];
					if (bases[mp2].protype == pclath) {
						dx = bases[mp].xcom - bases[mp2].xcom;
						dy = bases[mp].ycom - bases[mp2].ycom;
						dz = bases[mp].zcom - bases[mp2].zcom;
						r2 = dx * dx + dy * dy + dz * dz;
						if (r2 < tol2) {
							flag = 1;
							k = s2;
							j = s1;
							cout << " Cancel clathrin association, overlapping centers! " << mp << ' ' << mp2 << endl;
						}
					}
				}
			}
		}
	}

	return flag;

}
