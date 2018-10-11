#include "vector_rot_calls.h"

void rotate(double *v, double *M, double *r) {
	r[0] = M[0] * v[0] + M[1] * v[1] + M[2] * v[2];
	r[1] = M[3] * v[0] + M[4] * v[1] + M[5] * v[2];
	r[2] = M[6] * v[0] + M[7] * v[1] + M[8] * v[2];

}
