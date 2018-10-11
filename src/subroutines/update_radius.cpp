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

void update_radius(int c1, Complex *ind_com, Fullmol *bases) {
	/*For a multi-protein complex, determine how far from the COM is the 
	 farthest protein plus it's radius to define the complex's radius.
	 */
	int i;
	int p1;
	int size = ind_com[c1].mysize;
	ind_com[c1].radR = 0;
	double dx, dy, dz, mx, my, mz;
	double r1, r2;
	for (i = 0; i < size; i++) {
		p1 = ind_com[c1].plist[i];

		/*use the distance from the complex com to each protein com to find the largest radius
		 in each direction*/
		dx = ind_com[c1].xcom - bases[p1].xcom;
		dy = ind_com[c1].ycom - bases[p1].ycom;
		dz = ind_com[c1].zcom - bases[p1].zcom;
		r2 = dx * dx + dy * dy + dz * dz;
		r1 = sqrt(r2);
		mx = r1 + bases[p1].radR;
		if (mx > ind_com[c1].radR)
			ind_com[c1].radR = mx;

	}

}
