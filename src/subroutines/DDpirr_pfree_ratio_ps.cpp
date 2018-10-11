/*
 2D related 
 */

#include <stdio.h>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_vector.h>
#include "2Drelated.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <cfloat>

using namespace std;

double DDpirr_pfree_ratio_ps(gsl_matrix* mpir, gsl_matrix* msur, gsl_matrix* mnorm, double r, double Dtot, double deltat, double r0, double ps_prev, double rtol, double bindrad) {

	double pfree, pnormval, pirrval;
	double temp = r * r0 / (2.0 * Dtot * deltat);
	double temp2 = 1 / (4.0 * M_PI * deltat * Dtot);
	double temp3 = (exp(temp - (r0 * r0 + r * r) / (4.0 * deltat * Dtot)));
	/*RstepSize is coupled to the DDmatrixcreate RstepSize, they must be the SAME definition, 
	 *so RstepSize should not be defined any other way!
	 */
	const double RstepSize = sqrt(Dtot * deltat) / 50; 

	pfree = temp2 * temp3 * gsl_sf_bessel_I0_scaled(temp);
	pnormval = pnorm(mnorm, RstepSize, r0, bindrad);
	pirrval = pirr(mpir, msur, RstepSize, r, r0, bindrad);

	double pfreeN = pfree / pnormval; //NORMALIZES TO ONE

	double ratio;
	if (abs(pirrval - pfreeN * ps_prev) < rtol)
		ratio = 1.0;
	else {
		ratio = pirrval / (pfreeN * ps_prev);
	}
	return ratio;
}
