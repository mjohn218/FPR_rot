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

void TBLsur(gsl_matrix*& massc, double bindrad, double Dtot, double kr, double deltat, double Rmax, double RstepSize) {
	/////////////////////massc: the matrix for association probabilities///////////////////
	int ctr = 0;
	double result;//, RstepSize = sqrt(Dtot * deltat) / 50;
	//	double Rmax = 3.0 * sqrt(4.0 * Dtot * deltat) + bindrad;
	char funcID[] = "tblsur";

	gsl_function F;
	F.function = &fsur;
	f_params params;
	params.a = bindrad;
	params.D = Dtot;
	params.k = kr;
	params.t = deltat;

	gsl_set_error_handler_off();
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1e6);

	for (double Rindex = params.a; Rindex <= Rmax + RstepSize; Rindex += RstepSize) {
		params.r0 = Rindex;
		gsl_matrix_set(massc, 0, ctr, Rindex);
		F.params = reinterpret_cast<void *>(&params);

		if (kr < 1.0 / 0.0) {
			result = integrator(F, params, w, Rindex, bindrad, Dtot, kr, deltat, funcID, fsur);
			gsl_matrix_set(massc, 1, ctr, result);

		} else {
			result = integrator(F, params, w, Rindex, bindrad, Dtot, kr, deltat, funcID, fsur);
			gsl_matrix_set(massc, 1, ctr, 1.0 - result);
		}
		ctr = ctr + 1;
	}
	gsl_integration_workspace_free(w);
	gsl_set_error_handler(NULL);

}
