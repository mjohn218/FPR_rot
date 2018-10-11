#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void write_crds(Fullmol *bases, int p1) {
  cout << bases[p1].xcom << ' ' << bases[p1].ycom << ' ' << bases[p1].zcom << endl;
	int i;
	for (i = 0; i < bases[p1].ninterface; i++)
		cout << bases[p1].x[i] << ' ' << bases[p1].y[i] << ' ' << bases[p1].z[i] << endl;

}
