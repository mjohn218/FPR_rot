#include "reactions.h"

void ap2_coupledrxn_sub(int mu, int *Nmycoup, int **mycoupled, Fullmol *bases, int **Rlist, int *i_home, int *p_home, int p1, int p2) {
	/*If the ongoing reaction is coupled to the change in reaction states of other sites,
	 update those reactant species to their new reactions.
	 These coupled reactions are assumed zeroth-order, A->B 
	 and for r1>0 (regular), the change from A->B REMOVES an
	 interface that was previously free-to-bind.
	 but in special cases (where r1 is read in as <0)
	 For r1<0, the bound state of a molecule changes, but no free sites are removed.
	 */

	int i;
	int i1, prod, iind;
	double r1;
	int ptype;
	int wprot;
	int f, g;
	int curr, next;
	int oprot;
	for (i = 0; i < Nmycoup[mu]; i++) {
		r1 = mycoupled[mu][i];
		if (r1 < 0) {

			prod = Rlist[mu][2];
			ptype = p_home[prod];
			if (ptype == bases[p1].protype) {
				wprot = p1;

			} else {
				wprot = p2;

			}
			curr = Rlist[int(abs(r1))][0];
			next = Rlist[int(abs(r1))][1];
			for (f = 0; f < bases[wprot].ninterface; f++) {
				if (bases[wprot].istatus[f] == curr) {
					oprot = bases[wprot].partner[f];
					bases[wprot].istatus[f] = next;
					f = bases[wprot].ninterface;
					cout << "changing on p1: " << wprot << " rxn: " << r1 << " iface: " << curr << " to next: " << next << " partner: " << oprot << endl;
					for (g = 0; g < bases[oprot].nbnd; g++) {
						if (bases[oprot].bndlist[g] == curr) {
							if (bases[oprot].partner[g] == wprot) {
								bases[oprot].bndlist[g] = next;
								g = bases[oprot].nbnd;
								cout << "on partner, " << oprot << " changing iface: " << curr << " to next: " << next << endl;
							}
						}
					}
					for (g = 0; g < bases[oprot].ninterface; g++) {
						if (bases[oprot].istatus[g] == curr) {
							if (bases[oprot].partner[g] == wprot) {
								bases[oprot].istatus[g] = next;
								g = bases[oprot].ninterface;
							}
						}
					}
					/*done checking partner*/

				}
			}

			for (f = 0; f < bases[wprot].nbnd; f++) {
				if (bases[wprot].bndlist[f] == curr) {
					bases[wprot].bndlist[f] = next;
					f = bases[wprot].nbnd;
				}
			}

		} else {

			/*just regular coupled reaction
			 a free interface will be removed if there is a coupled reaction*/

		  i1 = Rlist[int(r1)][0];
		  prod = Rlist[int(r1)][1];
			iind = i_home[i1];
			ptype = p_home[i1];
			if (ptype == bases[p1].protype)
				wprot = p1;
			else
				wprot = p2;
			//remove this interface i1 from the freelist
			for (f = 0; f < bases[wprot].nfree; f++) {
				if (bases[wprot].freelist[f] == i1) {
					bases[wprot].freelist[f] = bases[wprot].freelist[bases[wprot].nfree - 1];
					f = bases[wprot].nfree;
				}
			}
			bases[wprot].nfree -= 1;
			bases[wprot].istatus[iind] = prod;

		}
	}
}
