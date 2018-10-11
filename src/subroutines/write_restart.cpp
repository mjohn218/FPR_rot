#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void write_restart(ofstream &outfile, Protein *wholep, Fullmol *bases, Parms &plist, int *Ncopy, int iter, Complex *ind_com, double deltat) {

	int i, j, p, n, nint,MAXPR=10;
	int ig;
	int f;
	int t = 0;
	
	char *names = new char[MAXPR];
	names[0] = 'A';
	names[1] = 'B';
	names[2] = 'C';
	names[3] = 'D';
	names[4] = 'E';
	names[5] = 'F';
	names[6] = 'G';
	names[7] = 'H';
	names[8] = 'I';
	names[9] = 'J';
	
	int Ntotalmol = plist.Ntotalmol;
	outfile << "total complexes: " << plist.ntotalcomplex << " iter: " << iter << endl;
	for (i = 0; i < Ntotalmol; i++) {
		outfile << "protein: " << i << " mycomplex: " << bases[i].mycomplex << endl;
		outfile << "Nfree: " << bases[i].nfree << " Freelist: ";
		for (f = 0; f < bases[i].nfree; f++)
			outfile << bases[i].freelist[f] << '\t';
		outfile << endl;
		outfile << "Nbound: " << bases[i].nbnd << " Boundlist: ";
		for (f = 0; f < bases[i].nbnd; f++)
			outfile << bases[i].bndlist[f] << '\t';
		outfile << endl;
		outfile << "Niface: " << bases[i].ninterface << " Istatus: ";
		for (j = 0; j < bases[i].ninterface; j++) {
			outfile << bases[i].istatus[j] << '\t';
		}
		outfile << endl;
		outfile << "Partners: ";
		for (j = 0; j < bases[i].ninterface; j++)
		  outfile << bases[i].partner[j] << '\t';
		outfile << endl;

	}
	/*From status of the proteins in the complex, should be able to reconstruct the
	 nature of each complex, its size, its prolist, and its diffusion and COM */

	outfile << plist.Natom << endl;
	outfile << "iter: " << iter << endl;

	for (i = 0; i < plist.Nprotypes; i++) {

		nint = wholep[i].ninterface;
		for (j = 0; j < Ncopy[i]; j++) {
			//outfile<<bases[t].xcom<<' '<<bases[t].ycom<<' '<<bases[t].zcom<<' '<<bases[t].mycomplex<<' '<<bases[t].mycomplexsize<<endl;
			outfile << std::fixed << std::setprecision(12) << names[i] << ' ' << bases[t].xcom << ' ' << bases[t].ycom << ' ' << bases[t].zcom << endl;
			for (n = 0; n < nint; n++) {
				outfile << std::fixed << std::setprecision(12) << names[i] << ' ' << bases[t].x[n] << ' ' << bases[t].y[n] << ' ' << bases[t].z[n] << endl;
			}
			t++;
		}
	}
	
	
}
