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

int break_sigma(int p1, int mu, int kind, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double *bindrad, int p2, int i1, int i2, int *p_home, int **myrxn) {

	/*Dissociate two bound interfaces, make sure all the proteins each is bound to
	 get put in the correct final complex.
	 proteins are already separated by sigma, so no further separation
	 */

	int c1 = bases[p1].mycomplex;
	int s1 = ind_com[c1].mysize;
	int prod = Rlist[mu][0];

	//change status of the interface
	int iind = ihome[i1];
	int iind2 = ihome[i2];
	int i, k;
	int mp, s2;
	int c2 = plist.ntotalcomplex;
	bases[p2].mycomplex = c2;

	/*assign each protein in original complex c1 to one of the two new complexes,
	 if the complex forms a loop, they will be put back together in c1, and the 
	 individual interfaces that dissociated freed.
	 */
	int cancel = 0;
	//	if (bases[p1].protype == bases[p2].protype && bases[p1].protype == plist.pclath)
	cancel = determine_which_complex_double(p1, p2, ind_com, c1, c2, bases, myrxn, Rlist, ihome, p_home, plist);
	//	else
	//cancel = determine_which_complex_merge(p1, p2, ind_com, c1, c2, bases, myrxn, Rlist, ihome, p_home, plist);
	
	if (cancel == 0) {

		/*continue on with the dissociation that creates two complexes*/

		bases[p1].istatus[iind] = i1;
		bases[p1].nbnd -= 1;
		bases[p1].bndlist[kind] = bases[p1].bndlist[bases[p1].nbnd];

		bases[p2].istatus[iind2] = i2; //used to be the product state
		bases[p2].nbnd -= 1;

		for (k = 0; k < bases[p2].nbnd + 1; k++) {
			if (bases[p2].bndlist[k] == prod) {
				bases[p2].bndlist[k] = bases[p2].bndlist[bases[p2].nbnd]; //replace in list with last
				k = bases[p2].nbnd + 1; //only replace one of the products, could have multiple 
			}
		}

		/*Add these protein into the bimolecular list*/
		bases[p1].nfree += 1;
		bases[p2].nfree += 1;
		bases[p1].freelist[bases[p1].nfree - 1] = i1;
		bases[p2].freelist[bases[p2].nfree - 1] = i2;

		/*create a new complex to hold partner and his partners*/
		plist.ntotalcomplex++;

		/*Update the NofEach Arrays for each new complex.*/
		update_NofEach(c1, plist, ind_com, bases);
		update_NofEach(c2, plist, ind_com, bases);
		/*Now we have the correct lists of proteins in complex 1 and complex2,
		 */
		update_diffusion(c1, ind_com, bases);
		update_diffusion(c2, ind_com, bases);

		update_one_com_only(c1, ind_com, bases);
		//update_com_only(c2, ind_com, bases);//will have to do this again after moving

		double stretch;
		double addx;
		double addy;
		double addz;
		double *chg1 = new double[3];
		double *chg2 = new double[3];

		double R = bindrad[mu];
		double Ry, R2, axe2;
		double sign = 1;
		
		double dx = bases[p1].x[iind] - bases[p2].x[iind2];
		
		
		double small = 1E-9;
		
		if (dx > 0)
		  dx = small;
		else
		  dx = -small; //so final separation is not binding radius with precision issues
		
		
		chg1[0] = dx; //Already separated by sigma, small adjustment for precision issues
		chg1[1] = 0;
		chg1[2] = 0;
		chg2[0] = 0;
		chg2[1] = 0;
		chg2[2] = 0;


		/*Update the positions of each protein. Then calculate the COMs of each
		 of the complexes. Then calculated the radius (uses the ind_com COM).
		 update rotational diffusion, translational diffusion was already updated.
		 */

		move_pro_coords(c1, ind_com, bases, chg1);
		move_pro_coords(c2, ind_com, bases, chg2);

		update_one_com_only(c1, ind_com, bases);
		update_one_com_only(c2, ind_com, bases);

		update_radius(c1, ind_com, bases);
		update_radius(c2, ind_com, bases);

		update_rot_diffusion(c1, ind_com, bases);
		update_rot_diffusion(c2, ind_com, bases);
		// cout <<"final crds: "<<bases[p1].xcom<<' '<<bases[p1].ycom<<' '<<bases[p1].zcom<<endl;
//     cout <<"final crds: "<<bases[p2].xcom<<' '<<bases[p2].ycom<<' '<<bases[p2].zcom<<endl;
//     cout <<"final cind_com: "<<ind_com[c1].xcom<<' '<<ind_com[c1].ycom<<' '<<ind_com[c1].zcom<<endl;

		delete[] chg1;
		delete[] chg2;

	} else {
		/*Reset all proteins back to complex c1, dissociation
		 will break the product state of the two proteins that dissociated but here they
		 are linked in a closed loop so it will not create a new complex.
		 positions don't change
		 */

		s1 = ind_com[c1].mysize;
		s2 = ind_com[c2].mysize;
		for (i = 0; i < s2; i++) {
		  mp = ind_com[c2].plist[i];
		  bases[mp].mycomplex = c1;
		  ind_com[c1].plist[s1] = mp;
		  s1++;
		}
		
		ind_com[c1].mysize = s1;

		/*Put interfaces that dissociated into free state*/
		bases[p1].istatus[iind] = i1;
		bases[p1].nbnd -= 1;
		bases[p1].bndlist[kind] = bases[p1].bndlist[bases[p1].nbnd];

		bases[p2].istatus[iind2] = i2; //used to be the product state
		bases[p2].nbnd -= 1;

		for (k = 0; k < bases[p2].nbnd + 1; k++) {
			if (bases[p2].bndlist[k] == prod) {
				bases[p2].bndlist[k] = bases[p2].bndlist[bases[p2].nbnd]; //replace in list with last
				k = bases[p2].nbnd + 1;
			}
		}

		/*Add these protein into the bimolecular list*/
		bases[p1].nfree += 1;
		bases[p2].nfree += 1;
		bases[p1].freelist[bases[p1].nfree - 1] = i1;
		bases[p2].freelist[bases[p2].nfree - 1] = i2;

	}

	return cancel;

}
