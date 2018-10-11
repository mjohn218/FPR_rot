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

int get_distance(Fullmol *bases, double **traj, int p1, int p2, double deltat, double bindrad, int *ncross, int **crosspart, int **crossint, int i1, int i2, int mu, int **cross_rxn, int *ihome, double iter, double Rmax, double &sep, double &R1, int *ncrosscom) {
	int i, j;

	double time = bases[p1].mytime;
	int i1ind = ihome[i1];
	int i2ind = ihome[i2];

	double dx = bases[p1].x[i1ind] - bases[p2].x[i2ind]; //these are original coordinates (pre movemovet)
	double dy = bases[p1].y[i1ind] - bases[p2].y[i2ind];
	double dz = bases[p1].z[i1ind] - bases[p2].z[i2ind];

	double R2 = dx * dx + dy * dy + dz * dz;
	R1 = sqrt(R2);
	sep = R1 - bindrad;

	int nc1 = ncross[p1];
	int nc2 = ncross[p2];
	double t1, t2;
	int flag = 0;
	/*Rmax should be the binding radius plus ~max diffusion distance, using 3*sqrt(6*Dtot*deltat)*/
	if (R1 < Rmax) {

		flag = 1;
		/*in this case we evaluate the probability of this reaction*/
		crosspart[p1][nc1] = p2;
		crosspart[p2][nc2] = p1;
		crossint[p1][nc1] = i1;
		crossint[p2][nc2] = i2;
		cross_rxn[p1][nc1] = mu;
		cross_rxn[p2][nc2] = mu;
		//		cout <<"cross: "<<p1<<' '<<p2<<" crossint: "<<i1ind<<" crossint2: "<<i2ind<<' '<<ncross[p1]<<' '<<ncross[p2]<<" cross interface: "<<i1<<' '<<i2<<" sep: "<<R1<<endl;
		ncross[p1]++;
		ncross[p2]++;
		ncrosscom[bases[p1].mycomplex]++;
		ncrosscom[bases[p2].mycomplex]++;

	} else {
		flag = 0;
	}
	return flag;
}
