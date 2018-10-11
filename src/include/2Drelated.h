#include <stdio.h>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_vector.h>

struct f_params { //Strucure for integrands
	double a; // sigma
	double D; // Diff. Coeff
	double k; // rate of assoc.
	double r0; // init. loc.
	double r; // init. rad.
	double t; // init. t.
	double rho; // current lipid density
};


void TBLpirr(gsl_matrix*& mpir, double bindrad, double Dtot, double kr, double deltat, double Rmax, double RstepSize);
void TBLsur(gsl_matrix*& msur, double bindrad, double Dtot, double kr, double deltat, double Rmax, double RstepSize);
void TBLnorm(gsl_matrix*& mnorm, double bindrad, double Dtot, double kr, double deltat, double Rmax, double RstepSize);


void DDmatrixcreate(gsl_matrix*& msur, gsl_matrix*& mnorm, gsl_matrix*& mpir, double bindrad, double Dtot, double kr, double deltat, double Rmax);
double DDpirr_pfree_ratio_ps(gsl_matrix* mpir, gsl_matrix* msur, gsl_matrix* mnorm, double r, double Dtot, double deltat, double r0, double ps_prev, double rtol, double bindrad);
int sizelookup(double bindrad, double Dtot, double deltat, double Rmax);
double pirr(gsl_matrix* mpir, gsl_matrix* msur, double RstepSize, double r, double r0, double a);
double DDpsur(gsl_matrix* msur, double Dtot, double deltat, double r0, double a);
double pnorm(gsl_matrix* mnorm, double RstepSize, double r0, double a);
double fsur(double x, void *p); //Gives us association probability not the survival
double fpir(double x, void *p);
double fpir_cum(double x, void *p);
double fnorm(double x, void *p);
double finfinity(double x, void *p);
double fassoc2dMF(double x, void *p);
double pirrfunc(double bindrad, double Dtot, double kr, double tval, double r0, double rcurr);
double assocfunct(double bindrad, double Dtot, double kr, double deltat, double r0);
double ktfunc(double bindrad, double Dtot, double kr, double deltat, double r0);
double callfassoc2dMF(double bindrad, double Dtot, double kr, double deltat, double rho);
int odefunc(double t, const double y[], double f[], void *params);
int jac(double t, const double y[], double *dfdy, double dfdt[], void *params);
double integrator(gsl_function F, f_params params, gsl_integration_workspace * w, double r0, double bindrad, double Dtot, double kr, double deltat, char* funcID, double (*f)(double, void*));
double ktintegrator(gsl_function F, f_params params, gsl_integration_workspace * w, double r0, double bindrad, double Dtot, double kr, double deltat, char* funcID, double (*f)(double, void*));


/*for exact GF sampling*/
void TBLpcume(gsl_matrix*& mpir, double bindrad, double Dtot, double kr, double deltat, double Rmax, double RstepSize, double Rlong);
void TBLpcume_allowerr(gsl_matrix*& mpir, int r0bins, int r1bins, double bindrad, double Dtot, double kr, double deltat, double Rmax, double RstepSize, double Rlong, gsl_matrix *pirr_mat);
void TBLpirr_longR(gsl_matrix*& mpir, double bindrad, double Dtot, double kr, double deltat, double Rmax, double RstepSize, double Rlong);
double integrate_numer(int r0ind, double bindrad, double rstar, gsl_matrix* pirr, double RstepSize, int Npt);
double integrate_trapz(int r0ind, double bindrad, double rstar, gsl_matrix* pirr, double RstepSize, int N);
