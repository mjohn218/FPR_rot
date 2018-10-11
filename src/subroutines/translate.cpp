#include "vector_rot_calls.h"

void translate(double *v, double *d, double *v2) {
	v2[0] = v[0] + d[0];
	v2[1] = v[1] + d[1];
	v2[2] = v[2] + d[2];

}
