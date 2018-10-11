#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "sphericalH.h"



using namespace std;

void gf_rot_1time(int Lmax, int thetabins, int azbins, double costheta, double phi, double time, double **prob, double Drot);

