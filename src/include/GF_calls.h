#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>

using namespace std;

double pirr_pfree_ratio_ps(double rcurr, double r0, double tcurr, double Dtot, double bindrad, double alpha, double ps_prev, double rtol);
double pirr_abs_pfree_ratio_ps(double rcurr, double r0, double tcurr, double Dtot, double bindrad, double ps_prev, double rtol);
double preflect_pfree_ratio(double rcurr, double r0, double tcurr, double Dtot, double bindrad, double alpharef, double rtol);

/*Use Faddeeva function to calculate exp()*erfc() when exp() argument is too large*/
double pirr_pfree_ratio_psF(double rcurr, double r0, double tcurr, double Dtot, double bindrad, double alpha, double ps_prev, double rtol);
double survive_irrF(double r0, double tcurr, double Dtot, double bindrad,double alpha, double cof);
double passocF(double r0, double tcurr, double Dtot, double bindrad,double alpha, double cof);

void pirrev_dist(int rbins, double r0, double time, double Dtot, double bindrad,double *pirrev, double alpha, double delr);
double survive_irr(double r0, double tcurr, double Dtot, double bindrad,double alpha, double cof);
double survive_irr_abs(double r0, double tcurr, double Dtot, double bindrad);
double passoc_absorbing(double r0, double tcurr, double Dtot, double bindrad);//ABSORBING BC in 3D

void pfree_dist(int rbins, double r0, double tcurr, double Dtot, double bindrad,double *pfree, double alpha, double delr, double passoc);
void pcorrection(int rbins, double delr, double deltat, double Dtot, double bindrad, double cof, double alpha_irr, double *prassoc, int r0bins, double delr0, double *pcorrect, double *pfreea);
void errormodel(double deltat, double Dtot, double bindrad, double kact, double &delr0, int r0bins, double *pcorrect, double &Rmax, double *pfreea);
double pirr_pfree_ratio(double rcurr, double r0, double deltat, double Dtot, double bindrad, double alpha_irr, double prevpassoc);
double pirr_pfree_diff(double rcurr, double r0, double deltat, double Dtot, double bindrad, double alpha_irr, double prevpassoc, double delr0);
double pirrev_value(double rcurr, double r0, double tcurr, double Dtot, double bindrad, double alpha);
double pfree_value(double rcurr, double r0, double tcurr, double Dtot, double bindrad,double alpha, double passoc);
double pfree_value_norm(double rcurr, double r0, double tcurr, double Dtot, double bindrad,double alpha);

double peval_cumulativeF(double r1, double r0, double tcurr, double Dtot, double bindrad, double alpha, double kact);
double peval_cumulative(double r1, double r0, double tcurr, double Dtot, double bindrad, double alpha, double kact);
double pirrev_valueF(double rcurr, double r0, double tcurr, double Dtot, double bindrad, double alpha);








