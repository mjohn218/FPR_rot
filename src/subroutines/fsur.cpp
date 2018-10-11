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

double fsur(double x, void *p) {

	double alp, bet, tet, P, T, h, f;
	f_params &params = *reinterpret_cast<f_params *>(p);

	if (params.k < 1.0 / 0.0) { //This gives the association probability
		h = (2.0 * M_PI * params.a * params.D);

		alp = h * x * j1(x * params.a) + params.k * j0(x * params.a);
		bet = h * x * y1(x * params.a) + params.k * y0(x * params.a);
		tet = alp * alp + bet * bet;
		P = j0(x * params.a) * y1(x * params.a) - j1(x * params.a) * y0(x * params.a);
		T = (j0(x * params.r0) * bet - y0(x * params.r0) * alp) / tet;

		f = T * P * (1 - exp(-params.D * params.t * x * x));
		f = f * params.a * params.k;
	} else { //absorbing boundary conditions... this gives the survival probability
		alp = j0(x * params.a);
		bet = y0(x * params.a);
		tet = alp * alp + bet * bet;
		P = j0(x * params.a) * y0(x * params.r0) - j0(x * params.r0) * y0(x * params.a);
		T = 1.0 / tet;

		f = (2.0 / M_PI) * T * P * exp(-params.D * params.t * x * x) / x;
	}

	return f;

}
