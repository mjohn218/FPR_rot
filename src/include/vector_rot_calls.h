#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>

using namespace std;

void dotproduct(double *v1, double *v2, double &theta);
void crossproduct(double *a, double *b, double *u);

/*rotation*/
void rotationEuler(double tx, double ty, double tz, double *M);
void rotationEulerX(double tx, double *M);
void rotationEulerY(double ty, double *M);
void rotationEulerZ(double tz, double *M);
void zrotationEuler(double tx, double ty, double *row);
void rotate(double *v, double *M, double *r);
void calc_Rmatrix(double *u, double theta, double *M);
void translate(double *v, double *d, double *v2);
void calcEuler(double &alpha, double &beta, double &gamma, double z1, double z2, double z3, double x3, double y3);
void rotationIntrinsic(double alpha, double beta, double gamma, double *M);
void rotationIntrinsicZYX(double alpha, double beta, double gamma, double *M);
void rotationEulerZXZ(double alpha, double beta, double gamma, double *M);
