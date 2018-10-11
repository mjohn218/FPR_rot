#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void write_protein_iface(ofstream &outfile, Parms &plist, Fullmol *bases, int *Ncopy, int it, Protein *wholep) {
	int i, j;
	int t = 0;
	outfile << plist.Natom << endl;
	outfile << "iter: " << it << endl;
	char zchar='Z';
	int n, nint;
	for (i = 0; i < plist.Nprotypes; i++) {

		nint = wholep[i].ninterface;
		for (j = 0; j < Ncopy[i]; j++) {
			//outfile<<bases[t].xcom<<' '<<bases[t].ycom<<' '<<bases[t].zcom<<' '<<bases[t].mycomplex<<' '<<bases[t].mycomplexsize<<endl;
		  outfile << std::fixed <<  zchar << ' ' << bases[t].xcom << ' ' << bases[t].ycom << ' ' << bases[t].zcom << endl;
		  //outfile << std::fixed << std::setprecision(12) << zchar << ' ' << bases[t].xcom << ' ' << bases[t].ycom << ' ' << bases[t].zcom << endl;
			for (n = 0; n < nint; n++) {
				outfile << std::fixed << zchar << ' ' << bases[t].x[n] << ' ' << bases[t].y[n] << ' ' << bases[t].z[n] << endl;
			}
			t++;
		}
	}

}
