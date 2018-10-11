#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "reactions.h"

#include "rand_gsl.h"
#include "Faddeeva.hh"
#include "utility_calls.h"
#include "vector_rot_calls.h"

void update_bound_proteins(int p1, int p2, int i1, int i2, Fullmol *bases,  Complex *ind_com, Parms &plist, int *ncrosscom, int prod, int iind, int iind2) 
{


  /*Combine the complexes into one complex by listing all c2 proteins into c1,
	  unless you are closing a loop, which could happen for clathrin cages.
	 Update the complex COM and diffusion constants
	 */
  int c1 = bases[p1].mycomplex;
  int s1 = ind_com[c1].mysize;
  int c2 = bases[p2].mycomplex;
  int s2 = ind_com[c2].mysize;
  int newsize = s1 + s2;
  
  int mp;
  double totalmassx = 0;
  double totalmassy = 0;
  double totalmassz = 0;
  ind_com[c1].xcom = 0;
  ind_com[c1].ycom = 0;
  ind_com[c1].zcom = 0;
  
  
    /*Update the status, free and bound lists of these two proteins*/
    //change status of the interface
    bases[p1].istatus[iind] = prod;
    bases[p2].istatus[iind2] = prod;
    //no longer free
    int i;
    for (i = 0; i < bases[p1].nfree; i++) {
      if (bases[p1].freelist[i] == i1) {
	bases[p1].freelist[i] = bases[p1].freelist[bases[p1].nfree - 1]; //put last reaction in place of this one
	i = bases[p1].nfree;
      }
    }
    bases[p1].nfree -= 1;
    for (i = 0; i < bases[p2].nfree; i++) {
      if (bases[p2].freelist[i] == i2) {
	bases[p2].freelist[i] = bases[p2].freelist[bases[p2].nfree - 1]; //put last reaction in place of this one
	i = bases[p2].nfree;
      }
    }
    bases[p2].nfree -= 1;
    
    /*add this as a possible dissociation reaction to both proteins
      so we know who they are bound to.
    */
    bases[p1].bndlist[bases[p1].nbnd] = prod; //put this reaction as last possible for uni
    bases[p1].nbnd += 1; //now a dissociation reaction is possible
    
    bases[p2].bndlist[bases[p2].nbnd] = prod;
    bases[p2].nbnd += 1; //now a dissociation reaction is possible
    
    bases[p1].partner[iind] = p2;
    bases[p2].partner[iind2] = p1;
    
    /*add c2's proteins to c1's list*/
    /*Get new COM of complex*/
    /*Copy final complex into the spot of now gone complex c2 */
    
    int j;
    
    if (c1 == c2) {
      cout << "CLOSING A LOOP! " << endl;
      plist.nloop++;
      update_complex_props(c1, plist, ind_com, bases);
    } else {
      /*Combine c1 and c2 into 1 complex, and replace c2 to shorten array.*/
      /*Update the NofEach species in the new complex*/
      for(int p=0;p<plist.Nprotypes;p++){
	int tmp=ind_com[c1].NofEach[p];
	ind_com[c1].NofEach[p]=tmp+ind_com[c2].NofEach[p];
      }
  
      ind_com[c1].mysize = newsize;
      /*Copy proteins from c2 into c1's list*/
      int t = s1;
      for (i = 0; i < s2; i++) {
	mp = ind_com[c2].plist[i];
	ind_com[c1].plist[t] = mp; //add these proteins to the first complex
	t++;
	bases[mp].mycomplex = c1;
      }
      update_complex_props(c1, plist, ind_com, bases);
      
      plist.ntotalcomplex -= 1;
      //cout <<"Before copying, complex c2 "<<c2<<" replacing with: "<<plist.ntotalcomplex<<endl;
      // write_crds_complex(c2, ind_com, bases);
      cout <<" before copying, coordinates to copy in: "<<plist.ntotalcomplex<<endl;
      //write_crds_complex(plist.ntotalcomplex, ind_com, bases);
      if (c2 != plist.ntotalcomplex) {
	
	/*otherwise, you are just deleting c2 entirely*/
	int tar = plist.ntotalcomplex;
	ncrosscom[c2] = ncrosscom[tar];
	ind_com[c2] = ind_com[tar];
	
	for (j = 0; j < ind_com[c2].mysize; j++) {
	  mp = ind_com[c2].plist[j];
	  bases[mp].mycomplex = c2;
	}
	cout <<"After copying, properties of c2, "<<c2<<"  should be same as  "<<plist.ntotalcomplex<< " above."<<endl;
	//write_crds_complex(c2, ind_com, bases);
      }
    }
    


}
