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
#include <time.h>
#include <sys/time.h>


using namespace std;

double integrator(gsl_function F, f_params params, gsl_integration_workspace * w, double r0, double bindrad, double Dtot, double kr, double deltat, char* funcID, double (*f)(double, void*)) {
  
	int it, status;
	int key=2;//21 point gauss kronod rules (max is 6. 
	double result, error, umax;
	const double xlow = 0, epsabs = 1e-7, epsrel = 1e-7;
	double publicationcrit = 1e-6;
	struct timeval clock;
	clock.tv_sec=0;
	clock.tv_usec=0;
	struct timeval duration;
	duration.tv_sec=0;
	duration.tv_usec=0;
	gettimeofday(&clock, NULL);
	/*For successful qagiu, takes ~0.01:0.3s, depending on R0. 
	  using qags even just twice takes ~4s-10s. qags seems necessary for large ka values. .*/
	status = gsl_integration_qagiu(&F, xlow, epsabs, epsrel, 1e6, w, &result, &error);
	

	it = 1;
	/*Try again with a weaker criterion*/
	if(status != GSL_SUCCESS){// && epsrel * (it + 1.0) < publicationcrit) {
	  //status = gsl_integration_qagiu(&F, xlow, epsabs * (it + 1.0), publicationcrit, 1e6, w, &result, &error);
	  status = gsl_integration_qagiu(&F, xlow, publicationcrit, publicationcrit, 1e6, w, &result, &error);
	  //	it += 1;
	  //	}
	}
	/*Now try with a different upper integration bound. Issues can arise due to oscillations of Bessel functions,
	  convergence can vary with choice of upper integration.
	 */
	if (status != GSL_SUCCESS) {
	  /*This value of umax is not optimized*/
	  umax = 10000.0;//sqrt(-log(DBL_MIN) / (Dtot * deltat));
	  while (abs((*f)(umax, &params)) > 1e-10) {
	    umax = umax * 1.2;
	  }
	  
	  cout << "______________" << endl;
		cout << funcID << " qagiu failed with status: " << status << endl;
		cout << "No solution found with rel/abs error smaller than: " << publicationcrit << endl;
		cout << "   @time: " << deltat << endl;
		cout << "   ka: " << kr << endl;
		cout << "   D: " << Dtot << endl;
		cout << "   r0: " << r0 << endl;
		cout << "Truncation will be performed on the semi-infinite domain..." << endl;

		it = 0;
		while (status != GSL_SUCCESS) {
		  /*Can we use qag, instead of qags?*/
		  status = gsl_integration_qag(&F, xlow, umax, epsabs, publicationcrit, 1e6, key, w, &result, &error);//epsabs was set to 0.0
			umax = umax * 0.9;
			it+=1;
		}

		cout << "New integration upper bound " << umax/0.9 << " found after " << it << " iterations." << endl;
		cout << "______________" << endl;
		struct timeval end_tv;
		gettimeofday(&end_tv, NULL);
		duration.tv_sec+=(end_tv.tv_sec-clock.tv_sec);
		duration.tv_usec+=(end_tv.tv_usec-clock.tv_usec);
		cout <<"Time of 2D table generation: "<<duration.tv_sec+1.0e-6*(double)duration.tv_usec<<endl;
		
	}
	return result;
}
