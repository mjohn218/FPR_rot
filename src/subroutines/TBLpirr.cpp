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

void TBLpirr(gsl_matrix*& mpir, double bindrad, double Dtot, double kr, double deltat, double Rmax, double RstepSize) {
	/////////////////////MPIRR////////////////////////
	int ctr1 = 0, ctr2 = 0;
	double result;
	//double RstepSize = sqrt(Dtot * deltat) / 50;
	//double Rmax = 3.0 * sqrt(4.0 * Dtot * deltat) + bindrad;
	char funcID[] = "tblpirr";

	gsl_function F;
	F.function = &fpir;
	f_params params;
	params.a = bindrad;
	params.D = Dtot;
	params.k = kr;
	params.t = deltat;

	gsl_set_error_handler_off();
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1e6);

	for (double Rindex = params.a; Rindex <= Rmax + RstepSize; Rindex += RstepSize) {

		params.r = Rindex;

		for (double R0index = params.a; R0index <= Rmax + RstepSize; R0index += RstepSize) {

			params.r0 = R0index;
			F.params = reinterpret_cast<void *>(&params);
			result = integrator(F, params, w, R0index, bindrad, Dtot, kr, deltat, funcID, fpir);
			gsl_matrix_set(mpir, ctr1, ctr2, result);
			ctr1 = ctr1 + 1;

		}
		ctr1 = 0;
		ctr2 = ctr2 + 1;
	}
	gsl_integration_workspace_free(w);
	gsl_set_error_handler(NULL);
}
