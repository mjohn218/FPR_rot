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

void update_complex_all(int Ncomplex, Complex *ind_com, Fullmol *bases) {
	/*Update Complex COM of all current complexes in
	 the system and also update their D and radius. 
	 */
	double totalmass = 0;

	int c1, s1;
	int mp;
	int i, j;
	for (c1 = 0; c1 < Ncomplex; c1++) {
		totalmass = 0;


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

		update_diffusion(c1, ind_com, bases);
		update_radius(c1, ind_com, bases);
	}

}
