#include "vector_rot_calls.h"

void calc_Rmatrix(double *u, double theta, double *M) {
	int i, j;
	double ux = u[0];
	double uy = u[1];
	double uz = u[2];
	double cthet = cos(theta);
	double sthet = sin(theta);
	M[0] = cthet + ux * ux * (1 - cthet);
	M[1] = ux * uy * (1 - cthet) - uz * sthet; //row 0, column 1, go across first!
	M[2] = ux * uz * (1 - cthet) + uy * sthet;
	M[3] = uy * ux * (1 - cthet) + uz * sthet;
	M[4] = cthet + uy * uy * (1 - cthet);
	M[5] = uy * uz * (1 - cthet) - ux * sthet;
	M[6] = uz * ux * (1 - cthet) - uy * sthet;
	M[7] = uz * uy * (1 - cthet) + ux * sthet;
	M[8] = cthet + uz * uz * (1 - cthet);

}
