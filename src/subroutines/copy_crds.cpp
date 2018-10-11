#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void copy_crds(int Nmol, Fullmol *bases, double *crd) {
	for (int i = 0; i < Nmol; i++) {
		crd[i * 3] = bases[i].xcom;
		crd[i * 3 + 1] = bases[i].ycom;
		crd[i * 3 + 2] = bases[i].zcom;
	}

}
