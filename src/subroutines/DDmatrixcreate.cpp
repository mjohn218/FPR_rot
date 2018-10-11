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

void DDmatrixcreate(gsl_matrix*& massc, gsl_matrix*& mnorm, gsl_matrix*& mpir, double bindrad, double Dtot, double kr, double deltat, double Rmax) {

	char fname[100];
	//	double Rmax=bindrad+3.0*sqrt(4.0*Dtot*deltat);
	double RstepSize=sqrt(Dtot*deltat)/50;

	TBLsur(massc, bindrad, Dtot, kr, deltat, Rmax, RstepSize);
	cout << "TBLassoc made successfully for Dtot:" << Dtot << ", Deltat:" << deltat << ", bindrad:" << bindrad << endl;

	TBLnorm(mnorm, bindrad, Dtot, kr, deltat, Rmax, RstepSize);
	cout << "TBLnorm made successfully for Dtot:" << Dtot << ", Deltat:" << deltat << ", bindrad:" << bindrad << endl;

	TBLpirr(mpir, bindrad, Dtot, kr, deltat, Rmax, RstepSize);
	cout << "TBLpirr made successfully for Dtot:" << Dtot << ", Deltat:" << deltat << ", bindrad:" << bindrad << endl;

	/*sprintf(fname, "MASSO_Dtot%3.3f_kr%3.3f_deltat%3.3f_bindrad%3.3f.dat", Dtot, kr, deltat, bindrad);
	FILE * f1 = fopen(fname, "w");
	gsl_matrix_fprintf(f1, massc, "%f");
	fclose(f1);

	sprintf(fname, "MPIR_Dtot%3.3f_kr%3.3f_deltat%3.3f_bindrad%3.3f.dat", Dtot, kr, deltat, bindrad);
	FILE * f2 = fopen(fname, "w");
	gsl_matrix_fprintf(f2, mpir, "%f");
	fclose(f2);

	sprintf(fname, "MNORM_Dtot%3.3f_kr%3.3f_deltat%3.3f_bindrad%3.3f.dat", Dtot, kr, deltat, bindrad);
	FILE * f3 = fopen(fname, "w");
	gsl_matrix_fprintf(f3, mnorm, "%f");
	fclose(f3);*/

}
