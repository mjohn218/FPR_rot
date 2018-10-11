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

double fnorm(double x, void *p) {

	double f;
	f_params &params = *reinterpret_cast<f_params *>(p);
	double temp = x * params.r0 / (2.0 * params.D * params.t);
	double temp2 = (x / (2.0 * params.D * params.t));
	double temp3 = (exp(temp - (params.r0 * params.r0 + x * x) / (4.0 * params.t * params.D)));
	f = temp2 * temp3 * gsl_sf_bessel_I0_scaled(temp);
	return f;
}
