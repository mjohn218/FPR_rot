#include "vector_rot_calls.h"

void crossproduct(double *a, double *b, double *u) {
	double u0 = a[1] * b[2] - a[2] * b[1];
	double u1 = a[2] * b[0] - a[0] * b[2];
	double u2 = a[0] * b[1] - a[1] * b[0];
	double r2 = u0 * u0 + u1 * u1 + u2 * u2;
	double r = sqrt(r2);
	if (r == 0)
		r = 1; //the vector will return 0 if they are same vector
	u[0] = u0 / r;
	u[1] = u1 / r;
	u[2] = u2 / r;

}
