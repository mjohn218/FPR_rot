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

double DDpsur(gsl_matrix* msur, double Dtot, double deltat, double r0, double a) {

	int index;
	double psurval;
	/*RstepSize is coupled to the DDmatrixcreate RstepSize, they must be the SAME definition,
 	 *so RstepSize should not be defined any other way!
	 */
	const double RstepSize = sqrt(Dtot * deltat) / 50;
	index = floor((r0 - a) / RstepSize);
	if (index < 0) {
		index = 0;
	}

	double val01 = gsl_matrix_get(msur, 1, index);
	double val02 = gsl_matrix_get(msur, 1, index + 1);
	double r01 = gsl_matrix_get(msur, 0, index);
	double r02 = gsl_matrix_get(msur, 0, index + 1);

	psurval = val01 * (r02 - r0) + val02 * (r0 - r01);
	psurval /= RstepSize;

	return psurval;

}
