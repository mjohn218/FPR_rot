#include "reactions.h"

void ap2_coupledrxn_add(int mu, int *Nmycoup, int **mycoupled, Fullmol *bases, int **Rlist, int *i_home, int *p_home, int p1, int p2) {
	/*If the ongoing reaction is coupled to the change in reaction states of other sites,
	 update those reactant species to their new reactions.
	 These coupled reactions are assumed zeroth-order, A->B 
	 and for r1>0 (regular), the change from A->B creates a new
	 interface that is free-to-bind.
	 but in special cases (where r1 is read in as <0)
	 For r1<0, the bound state of a molecule changes, but a new free site is not created.
	 */
	int i;
	int  i1, prod, iind;
	double r1;
	int ptype;
	int wprot;
	int f, g;
	int curr, next;
	int oprot;
	for (i = 0; i < Nmycoup[mu]; i++) {

		r1 = mycoupled[mu][i];

		if (r1 < 0) {
			/*Make pip2'2 bnd site 38 turn into a new site that doesn't dissociate*/

			prod = Rlist[mu][1];
			ptype = p_home[prod];
			if (ptype == bases[p1].protype) {
				wprot = p1;
				/*The other protein is p1'2 binding partner at the site below */
			} else {
				wprot = p2;

			}

			curr = Rlist[int(abs(r1))][0];
			next = Rlist[int(abs(r1))][1];
			for (f = 0; f < bases[wprot].ninterface; f++) {
				if (bases[wprot].istatus[f] == curr) {
					oprot = bases[wprot].partner[f]; //not the same as oprot above! 
					bases[wprot].istatus[f] = next;
					f = bases[wprot].ninterface;
					cout << "changing on p1: " << wprot << " rxn: " << r1 << " iface: " << curr << " to next: " << next << " partner: " << oprot << endl;
					/*Also need to change the bound state of the partner protein*/
					for (g = 0; g < bases[oprot].ninterface; g++) {
						if (bases[oprot].istatus[g] == curr) {
							//if(bases[oprot].partner[g]==wprot){
							/*Each protein should have only one copy of each bound state.*/
							bases[oprot].istatus[g] = next;
							g = bases[oprot].ninterface;
							//}
						}
					}
					for (g = 0; g < bases[oprot].nbnd; g++) {
						if (bases[oprot].bndlist[g] == curr) {
							//if(bases[oprot].partner[g]==wprot){THIS IS WRONG< partner[] is indexed by ninterface, not nbnd
							bases[oprot].bndlist[g] = next;
							g = bases[oprot].ninterface;
							cout << "on partner, " << oprot << " changing iface: " << curr << " to next: " << next << endl;
							//}
						}
					}

					/*done changing status of your bound partner*/

				}
			}

			for (f = 0; f < bases[wprot].nbnd; f++) {
				if (bases[wprot].bndlist[f] == curr) {
					bases[wprot].bndlist[f] = next; //bases[wprot].bndlist[bases[wprot].nbnd-1];
					f = bases[wprot].nbnd;
				}
			}

			/*For this reaction, no new interfaces opened up to bind*/
		} else {
			/*just regular coupled
			 there will be a new protein opened up as free to bind*/
		  i1 = Rlist[int(r1)][0];
		  prod = Rlist[int(r1)][1]; //free interface
			iind = i_home[prod];
			ptype = p_home[prod];

			if (ptype == bases[p1].protype)
				wprot = p1;
			else
				wprot = p2;

			bases[wprot].istatus[iind] = prod;
			bases[wprot].freelist[bases[wprot].nfree] = prod;
			bases[wprot].nfree += 1;
			//cout <<"Added free interface on protein: "<<wprot<<" interface: "<<prod<<" index: "<<iind<<endl;
		}
	}
}
