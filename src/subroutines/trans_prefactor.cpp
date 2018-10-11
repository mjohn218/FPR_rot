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

double trans_prefactor(double T, double nu, double scale, double a, double Dtot) {
	/*Dtrans=kbT/(6pinu a)
	 kb=1.3806488E-23m^2kg/(s^2 K)
	 T=293K
	 nu=0.001 kg/(ms) (at 293K) 1cP
	 --at 310K, nu is 0.000685 (or 0.685cP)
	 a is in m
	 Dtrans will be in m^2/s
	 so convert to read in a in nm, and Dtras return in nm^2/us
	 */
	double kb = 1.3806488E-23;
	double pre = kb * T / (6.0 * M_PI * nu) * 1E27 / 1E6; //1E27 is nm^3 and 1E6 is us
	/*Use scale factor if Einstein-Stokes Drot is 
	 changed to bead model Drot*/
	/*Use clathrin to scale*/
	//double a=bases[p1].radR;
	double ES = pre / a;
	cout << "Einstein-Stokes Dt: " << ES << " Input file Dt: " << Dtot << endl;
	if (scale > 0) {
	  scale = Dtot / ES; //correct for non spherical
	  pre*=scale;
	}

	return pre;

}
