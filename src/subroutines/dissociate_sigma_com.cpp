#include "reactions.h"
#include "rand_gsl.h"


void dissociate_sigma_com(int Ntotalmol, Fullmol *bases, Complex *ind_com, Parms &plist, int *p_home, int *i_home, int **Rlist, int *Nmyrxn, int **myrxn, double *bindrad, int *Ncoup, int **mycoupled, double irandmax, std::vector<int> &disslist, int *ncross, double *kr, int it, int *ncrosscom, int *movestat, int *nBoundPairs)
{
  
  int i, j;
  int dub;
  int n;
  int icom;
  int i1, i2, iind, iind2;
  int mu, ppart;
  double rate;
  double prob;
  double rnum, rnum2;
  int cancel;
  double us_to_s=1E-6;
  double deltat=plist.deltat;
  disslist.clear();
  for (i = 0; i < Ntotalmol; i++) {
    /*Test this proteins unimolecular interactions*/
    dub = 0;
    for (j = 0; j < bases[i].nbnd; j++) {
      icom = bases[i].bndlist[j]; //index of the bound species
      //cout <<"test dissociation : "<<i<<" Icom: "<<icom<<" iter: "<<it<<endl;
      for (n = 0; n < Nmyrxn[icom]; n++) {
	mu = myrxn[icom][n]; //reaction number
	
	i1 = Rlist[mu][1]; //determine interface for this reaction
	i2 = Rlist[mu][2];
	iind = i_home[i1];
	iind2=i_home[i2];
	if (p_home[i1] != p_home[i2]) {
	  if (p_home[i1] != bases[i].protype) {
	    i1 = Rlist[mu][2];//swap order of reactants
	    i2 = Rlist[mu][1];
	    iind = i_home[i1];
	    iind2 = i_home[i2];
	  }
	} else {
	  /*If two clathrins are dissociating, figure out which interface is on which protein.*/
	  if (i1 != i2) {
	    if (bases[i].istatus[iind] == icom) {
	      if (bases[i].istatus[iind2] == icom) {
		if (dub == 1) {
		  /*check other interface as well */
		  i1 = Rlist[mu][2];
		  i2 = Rlist[mu][1];
		  iind = i_home[i1];
		}
		dub++;
	      }
	    } else {
	      /*i1 is other interface */
	      i1 = Rlist[mu][2];
	      i2 = Rlist[mu][1];
	      iind = i_home[i1];
	      
	    }
	  } //if they are the same, no need to switch anything
	}
	
	rate = kr[mu];//kb
	ppart = bases[i].partner[iind];
	
	
	if (i < ppart) {
	  /*then this is the first protein in the reaction, so try it.
	    This if statement ensures we do not try to dissociate the same complex twice
	  */
	  
	  prob = 1 - exp(-rate * deltat * us_to_s);
	  rnum = 1.0 * rand_gsl();
	  
	  //cout <<"attempt dissociate: "<<i<<' '<<ppart<<" prob: "<<prob<<" rate: "<<rate<<" rnum: "<<rnum<<endl;
	  if (prob > rnum) {
	    rnum2 = rnum + rand_gsl() * irandmax; //to get higher resolution random numbers

	    if (prob > rnum2) {
	      cout << "Dissociate at iter: " << it << " protein: " << i << " partner: " << ppart << " randomnum: " << rnum2 << endl;
	      //								assemblyfile<<-1<<'\t'<<it*deltat<<'\t'<<i<<'\t'<<ppart<< endl;
	      /*Perform this dissociation reaction.
		Sometimes it is a bond broken, not a full dissociation to two complexes if
		the two interfaces are part of the same complex
	      */
	      
	      cancel = break_sigma(i, mu, j, bases, Rlist, i_home, ind_com, plist, bindrad, ppart, i1, i2, p_home, myrxn);
	      ap2_coupledrxn_sub(mu, Ncoup, mycoupled, bases, Rlist, i_home, p_home, i, ppart);
	      
	      //if(cancel==0)
	      update_Nboundpairs(bases[i].protype, bases[ppart].protype, -1, plist, nBoundPairs);
	      
	      // if (bases[i].protype == plist.pclath && bases[ppart].protype == plist.pclath){
	      // 	newlegcount=newlegcount-1; //updating the value of C-C/total legs
	      // 	cout <<"Dissociated 2 clathrins "<<endl;
	      // }else if (bases[i].protype == plist.pclath || bases[ppart].protype == plist.pclath) {
	      // 	newap2count=newap2count+1; /// updating the value of AP2-CLAT interaction. We are interested AP2 bounded/C-C
	      // 	cout <<"Dissociated AP2-clathrin " <<endl;
	      // }
	      cout << "Coords of p1 com: " << bases[i].xcom << "	" << bases[i].ycom << "	" << bases[i].zcom << endl;
	      cout << "Coords of p2 com: " << bases[ppart].xcom << "	" << bases[ppart].ycom << "	" << bases[ppart].zcom << endl;
	      
	      j--; //replaced this reaction, so stay on this one
	      disslist.push_back(i);
	      disslist.push_back(ppart);
	      
	      
	      /*If dissociated products are removed from overlap lists, use ncross=-1.
		If they remain in list to avoid overlap, use movestat=2 and also
		ensure that they are not allowed to diffuse again, by, for example,
		temporarily setting D=0.
	      */
	      //ncross[i] = -1;
	      //ncross[ppart] = -1;
	      //ncrosscom[bases[i].mycomplex]=-1;
	      //ncrosscom[bases[ppart].mycomplex]=-1;

	      
	      reflect_complex_rad_rot(i, bases, ind_com, plist);
	      reflect_complex_rad_rot(ppart, bases, ind_com, plist);
	      
	      //movestat[i]=2;
	      //movestat[ppart]=2;
	      
	      /*Need all proteins in each complex to be movestat=2, otherwise they will re-select their trajectories!
		With movestat=2, they will not be able to associate in this time-step. 
	       */
	      set_movestat_zero(i, bases, ind_com, movestat);
	      set_movestat_zero(ppart, bases, ind_com, movestat);
	      
	    }
	  }
	} //only try each pair dissociating once
      }
    } //finished trying out unimolecular reactions for this proteins
  }
  
}
