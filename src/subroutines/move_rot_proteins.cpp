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

void move_rot_proteins(int p1, Fullmol *bases, Complex *ind_com, double *traj, int *movestat, double *trajR, double *M) {
	int i, j;
	int k = bases[p1].mycomplex;
	int s1 = ind_com[k].mysize;
	double t1 = bases[p1].mytime;
	double *v = new double[3];
	double *v2 = new double[3];

	double dx = traj[0];
	double dy = traj[1];
	double dz = traj[2];
	rotationEuler(trajR[0], trajR[1], trajR[2], M);

	trajR[0] = 0;
	trajR[1] = 0;
	trajR[2] = 0;

	traj[0] = 0;
	traj[1] = 0;
	traj[2] = 0;
	
	double x0 = ind_com[k].xcom;
	double y0 = ind_com[k].ycom;
	double z0 = ind_com[k].zcom;

	ind_com[k].xcom += dx;
	ind_com[k].ycom += dy;
	ind_com[k].zcom += dz;

	//update protein COM
	int mp;
	for (i = 0; i < s1; i++) {
		mp = ind_com[k].plist[i];
		/*We've moved them, don't move them again*/

		movestat[mp] = 2;
		v[0] = bases[mp].xcom - x0;
		v[1] = bases[mp].ycom - y0;
		v[2] = bases[mp].zcom - z0;
		rotate(v, M, v2);
		/*first would make xcom=x0+v2, then would also add dx */
		bases[mp].xcom = x0 + dx + v2[0];
		bases[mp].ycom = y0 + dy + v2[1];
		bases[mp].zcom = z0 + dz + v2[2];

		//update interface coords
		for (j = 0; j < bases[mp].ninterface; j++) {
			v[0] = bases[mp].x[j] - x0;
			v[1] = bases[mp].y[j] - y0;
			v[2] = bases[mp].z[j] - z0;
			rotate(v, M, v2);
			/*first would make xcom=x0+v2, then would also add dx */
			bases[mp].x[j] = x0 + dx + v2[0];
			bases[mp].y[j] = y0 + dy + v2[1];
			bases[mp].z[j] = z0 + dz + v2[2];

		}
	}
	delete[] v;
	delete[] v2;
}
