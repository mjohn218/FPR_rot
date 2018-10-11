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

void TBLnorm(gsl_matrix*& mnorm, double bindrad, double Dtot, double kr, double deltat, double Rmax, double RstepSize) {
	////////////////////MNORM//////////////////////
	double result, error;
	const double xlowb = bindrad;
	int ctr = 0;
	const double epsabs = 1e-6;
	const double epsrel = 1e-6;
	//	double RstepSize = sqrt(Dtot * deltat) / 50;
	//double Rmax = 3.0 * sqrt(4.0 * Dtot * deltat) + bindrad;

	gsl_function F;
	F.function = &fnorm;
	f_params params;
	params.a = bindrad;
	params.D = Dtot;
	params.k = kr;
	params.t = deltat;

	for (double Rindex = params.a; Rindex <= Rmax + RstepSize; Rindex += RstepSize) {
		gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000000);
		params.r0 = Rindex;
		gsl_matrix_set(mnorm, 0, ctr, Rindex);
		F.params = reinterpret_cast<void *>(&params);
		gsl_integration_qagiu(&F, xlowb, epsabs, epsrel, 10000000, w, &result, &error);
		gsl_matrix_set(mnorm, 1, ctr, result);
		gsl_integration_workspace_free(w);
		ctr = ctr + 1;
	}
}
