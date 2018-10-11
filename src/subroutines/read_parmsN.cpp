#include "reactions.h"
#include "utility_calls.h"
#include <string>
#include <vector>

void read_parmsN(ifstream &parmfile, Parms &plist,std::vector<std::string> &infofilenames) {

	string tempstring[] = {"empty"};

	parmfile >> plist.Nprotypes;
	parmfile.ignore(400, '\n');
	parmfile >> plist.Nifaces;
	parmfile.ignore(400, '\n');
	parmfile >> plist.Nit;
	parmfile.ignore(400, '\n');
	parmfile >> plist.deltat;
	parmfile.ignore(400, '\n');
	parmfile >> plist.Nrxn;
	parmfile.ignore(400, '\n');
	parmfile >> plist.Nspecies;
	parmfile.ignore(400, '\n');
	parmfile >> plist.mass;
	parmfile.ignore(400, '\n');
	parmfile >> plist.xboxl;
	parmfile.ignore(400, '\n');
	parmfile >> plist.yboxl;
	parmfile.ignore(400, '\n');
	parmfile >> plist.zboxl;
	parmfile.ignore(400, '\n');
	parmfile >> plist.statwrite;
	parmfile.ignore(400, '\n');
	parmfile >> plist.configwrite;
	parmfile.ignore(400, '\n');
	parmfile >> plist.grwrite;
	parmfile.ignore(400, '\n');
	parmfile >> plist.restart;
	parmfile.ignore(400, '\n');
	parmfile >> plist.pclath;
	parmfile.ignore(400, '\n');

	parmfile >> tempstring[0];
	parmfile.ignore(400, '\n');
	for (int i=0;i<plist.Nprotypes;i++){
		parmfile >> tempstring[0];
		infofilenames.push_back(tempstring[0]);
	}
}
