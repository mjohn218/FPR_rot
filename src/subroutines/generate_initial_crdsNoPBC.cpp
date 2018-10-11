#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstring>
#include <sys/time.h>
#include "reactions.h"

#include "rand_gsl.h"
#include "Faddeeva.hh"
#include "utility_calls.h"
#include "vector_rot_calls.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

void generate_initial_crdsNoPBC(Parms plist, Fullmol *bases, int *Ncopy, Complex *ind_com, double *bindrad, Protein *wholep, int **Rlist,int *rxtype, int *p_home, int howmanylipids, double ***coordcont) {

	int i1,i2,i,j,jj,ii,k,MAXPR = 10,p,t = 0,mu=0, iface;
	int bflag = 1, Maxit = 50, bit = 0, noverlap = 0,rind, nint=0;
	int theyinteract = 0;
	double pbindr2, stretch, d2, dx, dy, dz,px,py,pz;
	double dx0,dx1,dy0,dy1,dz0,dz1,d0arm,d1arm;
	double d20,d21,d20arm,d21arm,dx0arm,dx1arm,dy0arm,dy1arm,dz0arm,dz1arm;
	char buffer[32], buffer2[32];
	char *let = new char[MAXPR];
	let[0] = 'A';
	let[1] = 'B';
	let[2] = 'C';
	let[3] = 'D';
	let[4] = 'E';
	let[5] = 'F';
	let[6] = 'G';
	let[7] = 'H';
	let[8] = 'I';
	let[9] = 'J';
	ofstream cfile("coordsALL.out");
	ofstream comfile("coordsCOM.out");

	// create random center of masses
	for (i = 0; i < plist.Ntotalmol; i++) {
		bases[i].xcom = plist.xboxl * rand_gsl() - plist.xboxl / 2.0;
		bases[i].ycom = plist.yboxl * rand_gsl() - plist.yboxl / 2.0;
		if (bases[i].Dz == 0) {
			bases[i].zcom = -plist.zboxl / 2.0;
		} else {
			bases[i].zcom = plist.zboxl * rand_gsl() - plist.zboxl / 2.0;
		}
	}

	while (bflag == 1 && bit < Maxit) {

		bit++;
		bflag = 0;
		noverlap = 0;

		for (i = 0; i < plist.Ntotalmol; i++) {
			for (j = howmanylipids; j < plist.Ntotalmol; j++) {
				if(i!=j){
					for(ii=0;ii<bases[i].ninterface;ii++){
						for(jj=0;jj<bases[j].ninterface;jj++){

							dx = bases[i].xcom + coordcont[bases[i].protype][ii+1][0] - (bases[j].xcom + coordcont[bases[j].protype][jj+1][0]);
							dy = bases[i].ycom + coordcont[bases[i].protype][ii+1][1] - (bases[j].ycom + coordcont[bases[j].protype][jj+1][1]);
							dz = bases[i].zcom + coordcont[bases[i].protype][ii+1][2] - (bases[j].zcom + coordcont[bases[j].protype][jj+1][2]);
							
							d2 = dx * dx + dy * dy + dz * dz;
							
							//do ii and jj interact? if so find mu
							theyinteract = 0;
							for (rind=0; rind<plist.Nrxn; rind++) {
								if(rxtype[rind] == 0){

									i1 = Rlist[rind][0];//interface1
									i2 = Rlist[rind][1];//interface2

									if ((bases[i].protype==p_home[i1] && bases[j].protype==p_home[i2]) || (bases[i].protype==p_home[i2] && bases[j].protype==p_home[i1])){
										theyinteract = 1;
										mu = rind;
										rind = plist.Nrxn;//exit loop
									}

								}
							}

							if(theyinteract==1){
								pbindr2 = bindrad[mu] * bindrad[mu];
							}
							else{
								pbindr2 = 0.0;
							}

							//make necessary measurements to check for overlap
							/*
						   *******
						   //This section here prevents the entire molecule from taking up space of other molecules, which can influence subsequent kinetics if the system is densely packed. It should not be enforced, because this type of volume exclusion (the whole molecule not justbinding sites) is not enforced as the simulation progresses!
							  dx0 = bases[i].xcom + coordcont[bases[i].protype][ii+1][0] - (bases[j].xcom);
							  dy0 = bases[i].ycom + coordcont[bases[i].protype][ii+1][1] - (bases[j].ycom);
							dz0 = bases[i].zcom + coordcont[bases[i].protype][ii+1][2] - (bases[j].zcom);

							d20 = dx0 * dx0 + dy0 * dy0 + dz0 * dz0;

							dx0arm = coordcont[bases[j].protype][jj+1][0];
							dy0arm = coordcont[bases[j].protype][jj+1][1];
							dz0arm = coordcont[bases[j].protype][jj+1][2];
							d20arm = dx0arm * dx0arm + dy0arm * dy0arm + dz0arm * dz0arm;

							dx1 = bases[i].xcom - (bases[j].xcom + coordcont[bases[j].protype][jj+1][0]);
							dy1 = bases[i].ycom - (bases[j].ycom + coordcont[bases[j].protype][jj+1][1]);
							dz1 = bases[i].zcom - (bases[j].zcom + coordcont[bases[j].protype][jj+1][2]);

							d21 = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;

							dx1arm = coordcont[bases[i].protype][ii+1][0];
							dy1arm = coordcont[bases[i].protype][ii+1][1];
							dz1arm = coordcont[bases[i].protype][ii+1][2];
							d21arm = dx1arm * dx1arm + dy1arm * dy1arm + dz1arm * dz1arm;

							if(d21<d21arm || d20<d20arm){

								bases[j].xcom = plist.xboxl * rand_gsl() - plist.xboxl / 2.0;
								bases[j].ycom = plist.yboxl * rand_gsl() - plist.yboxl / 2.0;
								if (bases[j].Dz == 0) {
									bases[j].zcom = -plist.zboxl / 2.0;
								} else {
									bases[j].zcom = plist.zboxl * rand_gsl() - plist.zboxl / 2.0;
								}
								noverlap++;
								bflag = 1;
								ii=bases[i].ninterface;
								jj=bases[j].ninterface;
//								cout<<i<<'\t'<<j<<endl;

}else*/
							if (d2 < pbindr2) {

								bases[j].xcom = plist.xboxl * rand_gsl() - plist.xboxl / 2.0;
								bases[j].ycom = plist.yboxl * rand_gsl() - plist.yboxl / 2.0;
								if (bases[j].Dz == 0) {
									bases[j].zcom = -plist.zboxl / 2.0;
								} else {
									bases[j].zcom = plist.zboxl * rand_gsl() - plist.zboxl / 2.0;
								}
								noverlap++;
								bflag = 1;
								ii=bases[i].ninterface;
								jj=bases[j].ninterface;
//								cout<<i<<'\t'<<j<<endl;

							}
						}
					}
				}

			}
		}
		cout << "Noverlap  " << noverlap << endl;
	}

	/*copy coords into interface pos and complex pos*/
	for (i = 0; i < plist.Ntotalmol; i++) {

		for (jj = 0; jj < bases[i].ninterface; jj++) {
			bases[i].x[jj] = bases[i].xcom+coordcont[bases[i].protype][jj+1][0];
			bases[i].y[jj] = bases[i].ycom+coordcont[bases[i].protype][jj+1][1];
			if (bases[i].Dz == 0) {
				bases[i].z[jj] = -plist.zboxl / 2.0 + coordcont[bases[i].protype][jj+1][2];
			} else {
				bases[i].z[jj] = bases[i].zcom+coordcont[bases[i].protype][jj+1][2];
			}
		}
		ind_com[i].xcom = bases[i].xcom;
		ind_com[i].ycom = bases[i].ycom;
		if (bases[i].Dz == 0) {
			ind_com[i].zcom = -plist.zboxl / 2.0;
		} else {
			ind_com[i].zcom = bases[i].zcom;
		}
		ind_com[i].plist[0] = i;
		ind_com[i].mysize = 1;
		bases[i].mycomplex = i;
		/*Add a species number tracker to each complex */
		ind_com[i].NofEach.reserve(plist.Nprotypes);
		for(int tmp=0;tmp<plist.Nprotypes;tmp++)
		  ind_com[i].NofEach.push_back(0);//fill array with zeros.
		ind_com[i].NofEach[bases[i].protype]=1;//has 1 copy of type bases[i].protype.
		//cout <<"i: "<<i<<" protype: "<<bases[i].protype<<" Complex N of Each Species, capacity "<<ind_com[i].NofEach.capacity()<<" Value: ";
		//for(int tmp=0;tmp<plist.Nprotypes;tmp++)
		//cout <<ind_com[i].NofEach[tmp]<<'\t';
		//cout <<endl;
		ind_com[i].Dx = bases[i].Dx;
		ind_com[i].Dy = bases[i].Dy;
		ind_com[i].Dz = bases[i].Dz;
		ind_com[i].Drx = bases[i].Drx;
		ind_com[i].Dry = bases[i].Dry;
		ind_com[i].Drz = bases[i].Drz;

	}

	t=0;
	for (p = 0; p < plist.Nprotypes; p++) {
		for (i = 0; i < Ncopy[p]; i++) {

			bases[t].radR = wholep[p].radR;
			ind_com[t].radR = wholep[p].radR;
			cfile << let[p] << " " << bases[t].xcom << ' ' << bases[t].ycom << ' ' << bases[t].zcom << endl;
			comfile << let[p] << " " << bases[t].xcom << ' ' << bases[t].ycom << ' ' << bases[t].zcom << endl;

			for (j = 0; j < wholep[p].ninterface; j++) {
				px = bases[t].xcom+coordcont[bases[t].protype][j+1][0];
				py = bases[t].ycom+coordcont[bases[t].protype][j+1][1];
				pz = bases[t].zcom+coordcont[bases[t].protype][j+1][2];
				sprintf(buffer2,"%d",j);
				cfile << let[p] <<buffer2<< " " << px << ' ' << py << ' ' << pz << endl;
			}
			t++;
		}
	}
	cout <<"Done generating initial coordinates "<<endl;
}
