#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "reactions.h"

#include "rand_gsl.h"
#include "Faddeeva.hh"
#include "utility_calls.h"
#include "vector_rot_calls.h"

void update_rot_diffusion(int c1, Complex *ind_com, Fullmol *bases) {
	int i, p1,j;
	int size = ind_com[c1].mysize;
	double Dxinv = 0;
	double Dyinv = 0;
	double Dzinv = 0;
	double inf = 1E300;
	p1 = ind_com[c1].plist[0];
	ind_com[c1].Drx = bases[p1].Drx;
	ind_com[c1].Dry = bases[p1].Dry;
	ind_com[c1].Drz = bases[p1].Drz;
	for (i = 0; i < size; i++) {

		p1 = ind_com[c1].plist[i];
		
		if (bases[p1].Drx != 0) {
		  Dxinv += 1.0 / pow(bases[p1].Drx, 1.0/3.0);
		} else {
			Dxinv = inf;
		}
		if (bases[p1].Dry != 0) {
		  Dyinv += 1.0 / pow(bases[p1].Dry, 1.0/3.0);
		} else {
			Dyinv = inf;
		}
		if (bases[p1].Drz != 0) {
		  Dzinv += 1.0 / pow(bases[p1].Drz, 1.0/3.0);
		} else {
			Dzinv = inf;
		}
	}
	ind_com[c1].Drx = 1.0 / (Dxinv*Dxinv*Dxinv);
	ind_com[c1].Dry = 1.0 /(Dyinv*Dyinv*Dyinv);
	ind_com[c1].Drz = 1.0 /(Dzinv*Dzinv*Dzinv);

	double tol=1E-50;
	if(ind_com[c1].Drx<tol)
	  ind_com[c1].Drx=0;
	if(ind_com[c1].Dry<tol)
	  ind_com[c1].Dry=0;
	if(ind_com[c1].Drz<tol)
	  ind_com[c1].Drz=0;
}

