#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>

using namespace std;


void cof_array(double **cofs, int Lmax, double **plmcof);
void cof_array_mpos(double **cofs, int Lmax);
double sqrt_fact_cof(int l, int m);
double fact_cof(int l, int m);
double lgamma_NM(double x);
void legendre(double **Plm, double x, int Lmax);
void Ylm(double **Plm, double **cof, double **plmcof, double **Yreal, double **Yimag, double phi, int Lmax);
void Ylm_mpos(double **Plm, double **cof, double **Yreal, double **Yimag, double phi, int Lmax);
