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

double pnorm(gsl_matrix* mnorm, double RstepSize, double r0, double a) {

	int index;
	double psurval;
	index = floor((r0 - a) / RstepSize);
	if (index < 0) {
		index = 0;
	}

	double val01 = gsl_matrix_get(mnorm, 1, index);
	double val02 = gsl_matrix_get(mnorm, 1, index + 1);
	double r01 = gsl_matrix_get(mnorm, 0, index);
	double r02 = gsl_matrix_get(mnorm, 0, index + 1);

	psurval = val01 * (r02 - r0) + val02 * (r0 - r01);
	psurval /= RstepSize;

	return psurval;

}
