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

int sizelookup(double bindrad, double Dtot, double deltat, double Rmax) {

	int ctr = 0;
	double stepSize = sqrt(Dtot * deltat) / 50;
	//	double Rmax = 3.0 * sqrt(4.0 * Dtot * deltat) + bindrad;
	for (double Rindex = bindrad; Rindex <= Rmax + stepSize; Rindex += stepSize) {
		ctr = ctr + 1;
	}
	return ctr;
}
