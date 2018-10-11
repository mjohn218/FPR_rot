#include "reactions.h"
#include "utility_calls.h"

int read_restartCELL(ifstream &restart, int Ntotalmol, Fullmol *bases, Parms &plist, Complex *ind_com, int *Ncopy, Protein *wholep) {
	int i, j, p, it, theit;
	int ig;
	int f;
	char str[300];
	cout << "RESTART FROM PRIOR CONFIG " << endl;
	//startfile.ignore(600,'\n');

	restart >> str >> str >> plist.ntotalcomplex >> str >> theit;
	for (i = 0; i < plist.ntotalcomplex; i++) {
		ind_com[i].mysize = 0;
	}
	int c1, n1, nfree, nbnd, niface;
	for (i = 0; i < Ntotalmol; i++) {
		restart >> str >> ig >> str >> c1;
		bases[i].mycomplex = c1;
		n1 = ind_com[c1].mysize;
		ind_com[c1].mysize++;
		ind_com[c1].plist[n1] = i;
		restart >> str >> nfree >> str;
		bases[i].nfree = nfree;
		for (j = 0; j < nfree; j++)
			restart >> bases[i].freelist[j];
		restart >> str >> nbnd >> str;
		bases[i].nbnd = nbnd;
		bases[i].npartner = bases[i].nbnd;
		for (j = 0; j < nbnd; j++)
			restart >> bases[i].bndlist[j];
		restart >> str >> niface >> str;

		for (j = 0; j < niface; j++)
			restart >> bases[i].istatus[j];
		restart >> str;
		for (j = 0; j < niface; j++)
			restart >> bases[i].partner[j];

	}

	cout << "total complexes: " << plist.ntotalcomplex << endl;
	for (i = 0; i < plist.ntotalcomplex; i++) {
		cout << "complex: " << i << " size: " << ind_com[i].mysize << endl;
	}

	int t = 0;
	restart >> str;
	restart >> str >> it;
	cout << str << endl;
	int n, nint;
	for (i = 0; i < plist.Nprotypes; i++) {
		nint = wholep[i].ninterface;
		for (j = 0; j < Ncopy[i]; j++) {
			restart >> str >> bases[t].xcom >> bases[t].ycom >> bases[t].zcom;
			cout << bases[t].xcom << "\t" << bases[t].ycom << "\t" << bases[t].zcom << endl;
			for (n = 0; n < nint; n++) {
				restart >> str >> bases[t].x[n] >> bases[t].y[n] >> bases[t].z[n];
				cout << bases[t].x[n] << "\t" << bases[t].y[n] << "\t" << bases[t].z[n] << endl;
			}
			t++;
		}
	}

	return theit;
}
