#include "reactions.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

void write_dcd(ofstream &outfile, Parms &plist, Fullmol *bases, int *Ncopy, int it, Protein *wholep, std::vector<std::string> &infofilenames) {
	int i, j;
	int t = 0;
//  outfile<<2*plist.Natomwrite-10<<endl;
//	outfile << plist.Natomwrite << endl;
//	outfile << "iter:" << it << endl;
	char aname;
	int n, nint, blx,bly,blz;
	int ind;
	double maxbondlengthsq = 3000.0;
	for (i = 0; i < plist.Nprotypes; i++) {

		nint = wholep[i].nint_write;
		for (j = 0; j < Ncopy[i]; j++) {
			//outfile<<bases[t].xcom<<' '<<bases[t].ycom<<' '<<bases[t].zcom<<' '<<bases[t].mycomplex<<' '<<bases[t].mycomplexsize<<endl;

//			outfile<<infofilenames[i].substr(0,3).c_str()<<j<<0<<' '<<bases[t].xcom<<' '<<bases[t].ycom<<' '<<bases[t].zcom<<endl;
			outfile<<bases[t].xcom<<' '<<bases[t].ycom<<' '<<bases[t].zcom<<endl;

			for (n = 0; n < nint; n++) {
				//ind=wholep[i].wrlist[n];
//				outfile << infofilenames[i].substr(0,3).c_str()<<j<<n+1<< ' ' << bases[t].x[n] << ' ' << bases[t].y[n] << ' ' << bases[t].z[n] << endl;
				if (wholep[i].Drx>0){
					blx = bases[t].x[n]-bases[t].xcom;
					bly = bases[t].y[n]-bases[t].ycom;
					blz = bases[t].z[n]-bases[t].zcom;
					if(blx*blx+bly*bly+blz*blz>maxbondlengthsq){
						outfile<<bases[t].xcom<<' '<<bases[t].ycom<<' '<<bases[t].zcom<<endl;
					}else{
						outfile << bases[t].x[n] << ' ' << bases[t].y[n] << ' ' << bases[t].z[n] << endl;
					}
				}else{
					outfile << bases[t].x[n] << ' ' << bases[t].y[n] << ' ' << bases[t].z[n] << endl;
				}
			}
			t++;
		}
	}

//  t = 0;
//  for(i=0;i<1;i++){
//
//    nint=wholep[i].nint_write;
//    for(j=0;j<Ncopy[i];j++){
//
//      for(n=0;n<nint;n++){
//    	  outfile<<'D'<<' '<<bases[t].x[n]<<' '<<bases[t].y[n]<<' '<<-200.0<<endl;
//      }
//      t++;
//    }
//  }

//  for(i=1;i<2;i++){
//
//    nint=wholep[i].nint_write;
//    for(j=0;j<Ncopy[i];j++){
//
//      for(n=0;n<nint;n++){
//    	  outfile<<'E'<<' '<<bases[t].x[n]<<' '<<bases[t].y[n]<<' '<<-200.0<<endl;
//      }
//      t++;
//    }
//  }

}
