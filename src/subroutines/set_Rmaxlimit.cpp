#include "reactions.h"
#include "utility_calls.h"
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

double set_Rmaxlimit(Parms plist, Protein *wholep, double *bindrad, int *i_home, int *p_home, int **Rlist, int *rxtype, double ***coordcont, double &Rmax_radius)
{

  /*For each reaction, need distance from the interface to the COM for both partners, plus the bindrad+sqrt(6*Dtot*deltat)
   */
  int j;
  double Rmaxlimit=0;
  double Rmaxtotal;
  double cf, r1C, r2C;
  Rmax_radius=0;
  cout << "N reactions: " << plist.Nrxn << endl;
  for (j = 0; j < plist.Nrxn; j++) {
    
    if (rxtype[j] == 0) {

      int i1=Rlist[j][0];
      int i2=Rlist[j][1];

      
      int pro1 = p_home[i1];
      int pro2 = p_home[i2];
      double sigma=bindrad[j];
      double Dtot=wholep[pro1].Dx+wholep[pro2].Dx;
      /*Add in contribution from Rotational diffusion.*/
      
      

      /*Now calculate distance from the interface to the protein COM.*/
      int iind1=i_home[i1];
      int iind2=i_home[i2];
      //      cout <<"In Set Rmaxlimit: Rxn: "<<j<<endl;
      cout <<"interface: "<<i1<<" index: "<<iind1<<" interface2: "<<i2<<" index: "<<iind2<<" pro1: "<<pro1<<" pro2: "<<pro2<<'\t';//" Dtot: "<<Dtot<<'\t';
      //cout <<"wholep[p1].Dx: "<<wholep[pro1].Dx<<" wholep[p2].Dx: "<<wholep[pro2].Dx<<endl;
      double dx=coordcont[pro1][iind1+1][0];//subtracted off COM already
      double dy=coordcont[pro1][iind1+1][1];//-coordcont[pro1][0][1];
      double dz=coordcont[pro1][iind1+1][2];//-coordcont[pro1][0][2];
      if(wholep[pro1].Dz==0){
	
	if(wholep[pro2].Dz==0)cout <<" IN RMAX LIMITE: 2D REACTION NUM: "<<j<<endl;
	double r1Csq=dx*dx+dy*dy;
	 r1C=sqrt(r1Csq);
        cf = cos(sqrt(2.0 * wholep[pro1].Drz * plist.deltat));
	double Dr1 = 2.0 * r1Csq * (1.0 - cf);
	Dtot+=Dr1/(4.0*plist.deltat);
      }else{
	double r1Csq=dx*dx+dy*dy+dz*dz;
	 r1C=sqrt(r1Csq);
	 cf = cos(sqrt(4.0 * wholep[pro1].Drz * plist.deltat));
	double Dr1 = 2.0 * r1Csq * (1.0 - cf);
	Dtot+=Dr1/(6.0*plist.deltat);
      }
      if(wholep[pro2].Dz==0){
	
	dx=coordcont[pro2][iind2+1][0];//-coordcont[pro2][0][0];
	dy=coordcont[pro2][iind2+1][1];//-coordcont[pro2][0][1];
	
	double r2Csq=dx*dx+dy*dy;
	 r2C=sqrt(r2Csq);
	cf = cos(sqrt(2.0 * wholep[pro2].Drz * plist.deltat));
	double Dr2= 2.0 * r2Csq * (1.0 - cf);
	
	Dtot+=Dr2/(4.0*plist.deltat);
      }else{
		
	dx=coordcont[pro2][iind2+1][0];//-coordcont[pro2][0][0];
	dy=coordcont[pro2][iind2+1][1];//-coordcont[pro2][0][1];
	dz=coordcont[pro2][iind2+1][2];//-coordcont[pro2][0][2];
	double r2Csq=dx*dx+dy*dy+dz*dz;
	 r2C=sqrt(r2Csq);
	cf = cos(sqrt(4.0 * wholep[pro2].Drz * plist.deltat));
	double Dr2= 2.0 * r2Csq * (1.0 - cf);
	
	Dtot+=Dr2/(6.0*plist.deltat);
      }
      
      double Rmax_diff=3.0*sqrt(6.0*Dtot*plist.deltat)+sigma;
      Rmaxtotal=Rmax_diff+r1C+r2C;
      cout <<"Rxn: "<<j<<" Rmax, for diffusion: "<<Rmax_diff<<" Deff, (including rotation): "<<Dtot<<" Interface1:COM dist: "<<r1C<<" Interface2:COM dist: "<<r2C<<" Total distaince: "<<Rmaxtotal<<endl;
      if(Rmaxtotal>Rmaxlimit){
	Rmaxlimit=Rmaxtotal;
	Rmax_radius=r1C+r2C;
      }
      
    }//only binding reactions

  }//all reactions
  cout <<"Rmaxlimit: "<<Rmaxlimit<<endl;
  return Rmaxlimit;
}
