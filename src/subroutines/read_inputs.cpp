#include "reactions.h"
#include "utility_calls.h"
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

void read_inputs(ifstream &rxnfile, Fullmol *bases, int *Ncopy, int *p_home, Protein *wholep, int *numpartner, int **Speclist, Parms &plist, int **Rlist, int *Ncoup, int **mycoupled, int *Nmyrxn, int **myrxn, double *bindrad, double *kr, int *cntrxn, int *freelist, int *bndlist, int *zlist, int *rxtype, int *ihome, std::vector<std::string> &infofilenames, double ***coordcont, int &howmanylipids) {

	/*Read in formatted list of reactions*/
	char rxntype;
	int i, pro1, pro2, npp, t=0, ct=0;
	int rxnnum, nint=0, iface;
	int p1, p2, p3;
	//char ctype[256];
	string ctype;
	
	int s;
	/*Only count reactions where you are a reactant!!!
	 as  a product, you just get created, you don't initiate any reactions 
	 */
	int j, k, flagexist=0;
	int nfree = 0;
	int nbnd = 0;
	int nmut = 0;
	double check2D=0.0;
	for (j = 0; j < plist.Nspecies; j++) {
	  Nmyrxn[j] = 0;
	}
	for (j = 0; j < plist.Nifaces; j++) {
		numpartner[j] = 0;
	}
	for (i = 0; i < plist.Nprotypes; i++) {
		wholep[i].npropart = 0;
	}

	//extract the structure from .info files
	for(i=0;i<plist.Nprotypes;i++){
		ifstream sfile(infofilenames[i].c_str()); // Use the constructor rather than `open`
		cout<<"Reading molecule info file "<<infofilenames[i]<<endl;
		if (sfile) // Verify that the file was opened successfully
		{
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');

			sfile >> wholep[i].Dx >> wholep[i].Dy >> wholep[i].Dz;
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');
			sfile >> wholep[i].Drx >> wholep[i].Dry >> wholep[i].Drz;
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');
			// sfile >> wholep[i].radx >> wholep[i].rady >> wholep[i].radz;
			// sfile.ignore(500, '\n');
			// sfile.ignore(500, '\n');
			// sfile.ignore(500, '\n');
			sfile >> nint; //number of protein i's interfaces
			wholep[i].ninterface = nint;
			wholep[i].nint_write = nint; //this can change below
			cout << "Protein: " << i << " Numinterfaces: " << wholep[i].ninterface << " nint value: "<<nint<<endl;
			cout << "Protein: " << i << " Dx: " << wholep[i].Dx << " Dz: " << wholep[i].Dz << endl;
			cout <<"CALCULATING THE RADIUS from coordinates, rather than from the input file!! "<<endl;
			//cout << "Protein: " << i << " radx: " << wholep[i].radx << " radz: " << wholep[i].radz << endl;
			for (j = 0; j < nint; j++) {
				sfile >> iface;
				cout << " index: " << iface << endl;
				wholep[i].valiface[j] = iface;
				p_home[iface] = i;
				ihome[iface] = j;
			}
			ct += nint;
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');
			sfile.ignore(500, '\n');
			if(wholep[i].Dz == 0){
				howmanylipids += Ncopy[i];
			}
			check2D+=wholep[i].Dz;
			//	initialize structure container to zero
			for(j=0;j<plist.Nifaces+1;j++){
				for (k=0;k<3;k++){
					coordcont[i][j][k] = 0.0;
				}
			}
			/*Read in all COM and interface coordinates*/
			for(j=0;j<wholep[i].ninterface+1;j++){
				sfile >> ctype;
				sfile >> coordcont[i][j][0] >> coordcont[i][j][1] >> coordcont[i][j][2];
				cout<< coordcont[i][j][0] << "\t"<< coordcont[i][j][1] << "\t"<< coordcont[i][j][2]<<endl;
				if(j>0){
					for(k=0;k<3;k++){
						coordcont[i][j][k] -= coordcont[i][0][k]; //make all interfaces with respect to the center of mass
					}
				}
			}
			/*Calculate the largest distance from the protein COM to interfaces*/
			wholep[i].radR=0;
			for(j=1;j<wholep[i].ninterface+1;j++){
			  double dx=coordcont[i][j][0];
			  double dy=coordcont[i][j][1];
			  double dz=coordcont[i][j][2];
			  
			  
			  double R1=sqrt(dx*dx+dy*dy+dz*dz);
			  cout <<" dx: "<<dx<<" dy: "<<dy<<" dz "<<dz<<" R1: "<<R1<<endl;
			  if(R1>wholep[i].radR){
			    wholep[i].radR=R1;
			  }
			}
			cout <<" PRotein: "<<i<<" Calculated radius: "<<wholep[i].radR<<endl;
			
			for (j = 0; j < Ncopy[i]; j++) {
			  
			  bases[t].ninterface = nint;//wholep[i].ninterface;
			  bases[t].protype = i;
			  bases[t].nfree = wholep[i].ninterface;

				bases[t].Dx = wholep[i].Dx;
				bases[t].Dy = wholep[i].Dy;
				bases[t].Dz = wholep[i].Dz;
				bases[t].Drx = wholep[i].Drx;
				bases[t].Dry = wholep[i].Dry;
				bases[t].Drz = wholep[i].Drz;

				bases[t].mass = wholep[i].radR; //plist.mass;
				if (wholep[i].radR == 0)
					bases[t].mass = 1;
				
				bases[t].nbnd = 0;
				bases[t].npartner = 0;
				bases[t].mycomplex = t;
				//cout <<"protein: "<<t<<" nint: "<<nint<<" Dx, Dz, Drx, Drz: "<<bases[t].Dx<<' '<<bases[t].Dz<<' '<<bases[t].Drx<<' '<<bases[t].Drz<<endl;
				t++;
			}

		}
		else
		{
			cerr << "File could not be opened!\n"; // Report error
			cerr << "Error code: " << strerror(errno); // Get some info as to why
			exit(1);
		}

	}

	rxnfile.ignore(600, '\n');
	cout << "N reactions: " << plist.Nrxn << endl;
	for (j = 0; j < plist.Nrxn; j++) {
		rxnfile >> rxnnum >> rxntype;
		if (rxntype == 'B') {
			rxtype[rxnnum] = 0;
			//binary reaction
			rxnfile >> Rlist[rxnnum][0] >> Rlist[rxnnum][1] >> Rlist[rxnnum][2];

			pro1 = p_home[Rlist[rxnnum][0]];
			pro2 = p_home[Rlist[rxnnum][1]];

			if (wholep[pro1].npropart<1){
				wholep[pro1].propart[wholep[pro1].npropart] = pro2;
				wholep[pro1].npropart++;
			}else{
				flagexist = 0;
				for(i=0;i<wholep[pro1].npropart;i++){
					if (wholep[pro1].propart[i]==pro2){
						flagexist = 1;
					}
				}
				if(flagexist==0){
					wholep[pro1].propart[wholep[pro1].npropart] = pro2;
					wholep[pro1].npropart++;
				}
			}
			Speclist[Rlist[rxnnum][0]][numpartner[Rlist[rxnnum][0]]] = Rlist[rxnnum][1];
			numpartner[Rlist[rxnnum][0]]++;

			if (Rlist[rxnnum][0]!=Rlist[rxnnum][1])	{
				if (wholep[pro2].npropart<1){
					wholep[pro2].propart[wholep[pro2].npropart] = pro1;
					wholep[pro2].npropart++;
				}else{
					flagexist = 0;
					for(i=0;i<wholep[pro2].npropart;i++){
						if (wholep[pro2].propart[i]==pro1){
							flagexist = 1;
						}
					}
					if(flagexist==0){
						wholep[pro2].propart[wholep[pro2].npropart] = pro1;
						wholep[pro2].npropart++;
					}
				}
				Speclist[Rlist[rxnnum][1]][numpartner[Rlist[rxnnum][1]]] = Rlist[rxnnum][0];
				numpartner[Rlist[rxnnum][1]]++;
			}

			rxnfile >> Ncoup[rxnnum];
			for (i = 0; i < Ncoup[rxnnum]; i++) {
				rxnfile >> mycoupled[rxnnum][i];
			}

			p1 = Rlist[rxnnum][0];
			p2 = Rlist[rxnnum][1];
			s = Nmyrxn[p1];
			myrxn[p1][s] = rxnnum;
			Nmyrxn[p1]++;
			freelist[nfree] = p1;
			nfree++;
			if (p1 != p2) {
				s = Nmyrxn[p2];
				myrxn[p2][s] = rxnnum;
				Nmyrxn[p2]++;
				freelist[nfree] = p2;
				nfree++;
			}

		} else if (rxntype == 'U') {
			//unimolecular reaction
			rxtype[rxnnum] = 1;
			rxnfile >> Rlist[rxnnum][0] >> Rlist[rxnnum][1] >> Rlist[rxnnum][2];
			rxnfile >> Ncoup[rxnnum];
			for (i = 0; i < Ncoup[rxnnum]; i++) {
				rxnfile >> mycoupled[rxnnum][i];
			}
//			/*For unbinding, also need to read in Kd in uM!!!!!!*/
//			rxnfile >> Kd[rxnnum];
			p1 = Rlist[rxnnum][0];
			s = Nmyrxn[p1];
			myrxn[p1][s] = rxnnum;
			Nmyrxn[p1]++;
			bndlist[nbnd] = p1;
			nbnd++;

		} else if (rxntype == 'Z') {
			//zeroth order
			rxtype[rxnnum] = 2;
			rxnfile >> Rlist[rxnnum][0] >> Rlist[rxnnum][1];
			rxnfile >> Ncoup[rxnnum];
			for (i = 0; i < Ncoup[rxnnum]; i++) {
				rxnfile >> mycoupled[rxnnum][i];
			}
			p1 = Rlist[rxnnum][0];
			zlist[nmut] = p1;
			nmut++;

		}
		/*Need reaction rates for the dissociation reactions*/
		rxnfile >> bindrad[rxnnum] >> kr[rxnnum];
//		cout << bindrad[rxnnum] << ' '<< kr[rxnnum]<< endl;
	}
	cntrxn[0] = nfree;
	cntrxn[1] = nbnd;
	cntrxn[2] = nmut;
	t=0;
	for (int p = 0; p < plist.Nprotypes; p++) {
		npp = wholep[p].npropart;

		for (i = 0; i < Ncopy[p]; i++) {
			bases[t].npropart = npp;
			for (j = 0; j < npp; j++)
				bases[t].propart[j] = wholep[p].propart[j];

			for (j = 0; j < wholep[p].ninterface; j++) {
				bases[t].istatus[j] = wholep[p].valiface[j];
				bases[t].freelist[j] = wholep[p].valiface[j];
				bases[t].partner[j]=-1;//initialize this variable, it is only set if the protein binds.
			}
			t++;
		}
	}

	for (i = 0; i < plist.Nprotypes; i++) {
		cout << "Protein: " << i << " Npropart: " << wholep[i].npropart << " Partners:" << endl;
		for (j = 0; j < wholep[i].npropart; j++) {
			cout << wholep[i].propart[j] << endl;
		}
	}

	for (i = 0; i < plist.Nifaces; i++) {
		cout << "Interface: " << i << " Nintpartner: " << numpartner[i] << " Partners:" << endl;
		for (j = 0; j < numpartner[i]; j++) {
			cout << Speclist[i][j] << endl;
		}

	}

	if (ct != plist.Nifaces) {
		cerr << "Number of interfaces assigned to proteins does not match network interface numbers!" << endl;
		exit(1);
	}

	if(check2D==0){
	  howmanylipids=0;
	  cout <<"---------------------------------"<<endl;
	  cout <<"SYSTEM IS IN 2D "<<endl;
	}
}
