/*
Generate a psf file for use with VMD
 */
#include "reactions.h"
#include "utility_calls.h"
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <stdio.h>

using namespace std;

void gen_psf_system(Parms &plist, Protein *wholep, int *Ncopy,std::vector<std::string> &infofilenames)
{
	//	output file info
	FILE * outfile;
	char fname[100];
	sprintf(fname, "system.psf");
	outfile = fopen (fname,"w");

	std::string atomnames = "OZNCSPF";

	int i, j, k, Ntotal=0, Nbond=0, t=1, zero=0, resnum=0, type=0, ncomm=2, n=0;
	double charge=0.5, mass=10.0;

	for(i=0;i<plist.Nprotypes;i++){
		Ntotal += (wholep[i].ninterface+1)*Ncopy[i];
		Nbond += wholep[i].ninterface*Ncopy[i];
	}

	//Write PSF Header
	fprintf(outfile, "PSF CMAP CHEQ\n");
	fprintf(outfile,"\n");
	fprintf(outfile,"%8d !NTITLE\n", ncomm);
	fprintf(outfile, "* PSF for whole System,  \n");
	fprintf(outfile,"* PSF comment total molecule %d\n", Ntotal);
	fprintf(outfile,"\n");
	fprintf(outfile,"%8d !NATOM\n",Ntotal);

	std::string segid = "CLA";

	for(i=0; i<plist.Nprotypes; i++){
		for(j=0; j<Ncopy[i]; j++){
			for(k=0; k<wholep[i].ninterface+1; k++){
				fprintf(outfile,"%8d %-4s %-4d %-4s %-4s %4d %10.6f    %10.4f  %10d\n", t, infofilenames[i].substr(0,3).c_str(), resnum, infofilenames[i].substr(0,3).c_str(), atomnames.substr(i,1).c_str(), type, charge, mass, zero);
				t++;
			}
			resnum++;
		}
	}

	//Write Bond Info
	t=1;
	fprintf(outfile,"\n");
	fprintf(outfile,"\n");
	fprintf(outfile,"%8d !NBOND: bonds\n", Nbond);

	for(i=0;i<plist.Nprotypes;i++){
		for(j=0; j<Ncopy[i]; j++){
			for(k=0;k<wholep[i].ninterface;k++){
				fprintf(outfile, "%8d%8d",t, t+k+1);
				n++;
				if(n%4==0)fprintf(outfile,"\n");
			}

			t+=wholep[i].ninterface+1;

		}
	}

	fclose(outfile);

}
