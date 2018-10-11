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

void update_diffusion(int c1, Complex *ind_com, Fullmol *bases) {

	/*Update the diffusion of a complex bases on sum of their radii,
	 defined via the Einstein-stokes equation, such that diffusions sum via their
	 inverses.
	 */
	int i, p1, j;
	int size = ind_com[c1].mysize;
	double Dxinv = 0;
	double Dyinv = 0;
	double Dzinv = 0;
	double inf = 1E300;
	p1 = ind_com[c1].plist[0];
	double currDx=ind_com[c1].Dx;
	double currDy=ind_com[c1].Dy;
	double currDz=ind_com[c1].Dz;

	ind_com[c1].Dx = bases[p1].Dx;
	ind_com[c1].Dy = bases[p1].Dy;
	ind_com[c1].Dz = bases[p1].Dz;
	for (i = 0; i < size; i++) {

		p1 = ind_com[c1].plist[i];
		//Dxinv+=1.0/bases[p1].Dx;
		//Dyinv+=1.0/bases[p1].Dy;

		if (bases[p1].Dx != 0) {
			Dxinv += 1.0 / bases[p1].Dx;
		} else {
			Dxinv = inf;
			//      ind_com[c1].Dx=0;
		}
		if (bases[p1].Dy != 0) {
			Dyinv += 1.0 / bases[p1].Dy;
		} else {
			Dyinv = inf;
			//ind_com[c1].Dy=0;
		}
		if (bases[p1].Dz != 0) {
			Dzinv += 1.0 / bases[p1].Dz;
		} else {
			Dzinv = inf;
			//ind_com[c1].Dz=0;
		}
	}
	ind_com[c1].Dx = 1.0 / Dxinv;
	ind_com[c1].Dy = 1.0 / Dyinv;
	ind_com[c1].Dz = 1.0 / Dzinv;
	double tol=1E-50;
	if(ind_com[c1].Dx<tol)
	  ind_com[c1].Dx=0;
	if(ind_com[c1].Dy<tol)
	  ind_com[c1].Dy=0;
	if(ind_com[c1].Dz<tol)
	  ind_com[c1].Dz=0;

	if(ind_com[c1].Dz<tol){
	  /*On the membrane. Only allow 2D diffusion at certain intervals, to avoid generating too many 2D Tables.*/
	  double dtmp;
	  if(ind_com[c1].Dx<0.0001)
	    dtmp=ind_com[c1].Dx*100000;
	  else if(ind_com[c1].Dx<0.001)
	    dtmp=ind_com[c1].Dx*10000;
	  else if(ind_com[c1].Dx<0.01)
	    dtmp=ind_com[c1].Dx*1000;
	  else if(ind_com[c1].Dx<0.1)
	    dtmp=ind_com[c1].Dx*100;
	  else
	    dtmp=ind_com[c1].Dx*10;
	  /*Keep only one sig fig for <1, 2 for 1<d<10, 3 for 10<d<100, etc*/
	  int d_ones=int(round(dtmp));
	  cout <<"FORCED D_2D to fewer sig-figs, starting value: "<<ind_com[c1].Dx<<" final value: ";
	  /*Now put back in correct size*/
	  if(ind_com[c1].Dx<0.0001)
	    ind_com[c1].Dx=d_ones*0.00001;
	  else if(ind_com[c1].Dx<0.001)
	    ind_com[c1].Dx=d_ones*0.0001;
	  else if(ind_com[c1].Dx<0.01)
	    ind_com[c1].Dx=d_ones*0.001;
	  else if(ind_com[c1].Dx<0.1)
	    ind_com[c1].Dx=d_ones*0.01;
	  else
	    ind_com[c1].Dx=d_ones*0.1;

	  cout <<ind_com[c1].Dx<<endl;
	  ind_com[c1].Dy=ind_com[c1].Dx;//set Dy equal to Dx.
	}
}
