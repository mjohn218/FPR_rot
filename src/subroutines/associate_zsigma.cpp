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

int associate_zsigma(int p1, int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom) 
{

	/*Move the interfaces that are reacting with one another to
	 the position geometrically between them, but move each
	 interface along this vector proportional to their diffusion constant
	 In this version, put the two interfaces into the same z-plane before calculating
	 the vector between them so they are not rotated around z.
	 So the molecules are adjacent to one another.
	 
	 No rotation of the complexes occurs, only translation.
	 */

	int prod = Rlist[mu][2];
	int iind = ihome[i1];
	int iind2 = ihome[i2];
	
	/*Store each protein's pre-association complex values*/
	int c1 = bases[p1].mycomplex;
	int c2 = bases[p2].mycomplex;

	int cancel;

	double *v0 = new double[3];
	double *v2 = new double[3];
	double *M = new double[9];
	double *u = new double[3];
	double *dtrans = new double[3];
	double *drev = new double[3];

	/*Now move protein*/
	double Dxtot = ind_com[c1].Dx + ind_com[c2].Dx;
	double Dytot = ind_com[c1].Dy + ind_com[c2].Dy;
	double Dztot = ind_com[c1].Dz + ind_com[c2].Dz;

	/*For proteins bound to membrane, Dztot could be zero*/
	double tol = 1E-16;
	if (Dztot < tol)
		Dztot = 1; //otherwise you divide by zero

	double dz = bases[p2].z[iind2] - bases[p1].z[iind];
	//dz -= plist.zboxl * round(dz / plist.zboxl);

	/*eliminate the differences in their z-position*/
	double delz1 =  dz * ind_com[c1].Dz / Dztot;// * (bindrad / dz - 1.0);
	double delz2 = -dz * ind_com[c2].Dz / Dztot;// * (bindrad / dz - 1.0);
	
	move_zcrds(bases, p1, p2, delz1, delz2, ind_com);

	/*The sign ensures the COM's get rotated
	 to the outside of interfaces final spot*/
	v0[0] = -bases[p2].x[iind2] + bases[p1].x[iind];
	
	v0[1] = -bases[p2].y[iind2] + bases[p1].y[iind];
	
	v0[2] = -bases[p2].z[iind2] + bases[p1].z[iind]; //This should be zero now!
	v0[2] = 0; //it is possible this is not zero if both associating molecules were bound to the membrane already
	
	 cout <<"After moving in z, current separtion: "<<endl;
 	cout <<" dx: "<<v0[0]<<' '<<v0[1]<< ' '<<v0[2]<<endl;
// 	cout <<"protein 1: "<<endl;
	//	write_crds(bases, p1);

// 	cout <<"protein 2: "<<endl;
 	//write_crds(bases, p2);

	double R2 = v0[0] * v0[0] + v0[1] * v0[1];
	double R1 = sqrt(R2);

	double theta;
	double pivx;
	double pivy;
	double pivz;
	double small = 1E-9;
	//distance to move to place interfaces, along vector v
	dtrans[0] = v0[0] * ind_com[c1].Dx / Dxtot * (bindrad / R1 - 1.0);
	dtrans[1] = v0[1] * ind_com[c1].Dy / Dytot * (bindrad / R1 - 1.0);
	dtrans[2] = v0[2] * ind_com[c1].Dz / Dztot * (bindrad / R1 - 1.0);

	drev[0] = -v0[0] * ind_com[c2].Dx / Dxtot * (bindrad / R1 - 1.0);
	drev[1] = -v0[1] * ind_com[c2].Dy / Dytot * (bindrad / R1 - 1.0);
	drev[2] = -v0[2] * ind_com[c2].Dz / Dztot * (bindrad / R1 - 1.0); //should be zero

	//cout <<"dtrans: "<<dtrans[0]<<' '<<dtrans[1]<<' '<<dtrans[2]<<endl;
	//cout <<"drev: "<<drev[0]<<' '<<drev[1]<<' '<<drev[2]<<endl;
	
	/*NOW WE ARE MOVING THE TWO PROTEINS TO A SEPARATION OF SIGMA*/
	//just translate if larger complexes are coming together
	translate_int(p1, c1, ind_com, bases, dtrans);

	/*Same for complex2*********************************/
	//just translate
	translate_int(p2, c2, ind_com, bases, drev);

	cout <<"Initial separation, no z :"<<R1<<'\t';
	// v0[0] = -bases[p2].x[iind2] + bases[p1].x[iind];
// 	v0[0] -= plist.xboxl * round(v0[0] / plist.xboxl);
// 	v0[1] = -bases[p2].y[iind2] + bases[p1].y[iind];
// 	v0[1] -= plist.yboxl * round(v0[1] / plist.yboxl);
// 	v0[2] = -bases[p2].z[iind2] + bases[p1].z[iind]; //This should be zero now!
// 	v0[2] = 0; //it is possible this is not zero if both associating molecules were bound to the membrane already
// 	R2 = v0[0] * v0[0] + v0[1] * v0[1];
// 	R1 = sqrt(R2);
	// cout <<"final sep: "<<R1<<endl;
	cout <<"final crds: "<<bases[p1].xcom<<' '<<bases[p1].ycom<<' '<<bases[p1].zcom<<endl;
	cout <<"final crds: "<<bases[p2].xcom<<' '<<bases[p2].ycom<<' '<<bases[p2].zcom<<endl;
	
	/*Determine if you crashed two clathrins together
	  unless they are in the same complex, in which case you are closing a loop.
	*/
	cancel=0;
	if (c1 != c2)
	  cancel = measure_overlap(c1, c2, ind_com, bases, plist.pclath);
	if(cancel==0)
	  cancel=check_box_span(p1, p2, bases, ind_com, plist);//if complex is too big to fit in the box, reverse the association!

	
	if(cancel==0)
	  update_bound_proteins( p1,  p2,  i1,  i2,  bases,  ind_com,  plist,  ncrosscom, prod,  iind, iind2);
	else{
	  
	  //  cout <<"ind_com crds: "<<ind_com[c1].xcom<<' '<<ind_com[c1].ycom<<' '<<ind_com[c1].zcom<<endl;
	  
	  /*CANCEL=1: IN THIS CASE< CLATHRINS CRASHED TOGETHER
	    un-bind them, allow to diffuse apart*/
	  cout << "Unbind clathrins, crashed together ! " << endl;
	  
	  dtrans[0] *= -1;
	  dtrans[1] *= -1;
	  dtrans[2] *= -1;
	  
	  drev[0] *= -1;
	  drev[1] *= -1;
	  drev[2] *= -1;
	  
	  /*Update the position of proteins in complex one*/
	  translate_int(p1, c1, ind_com, bases, dtrans);
	  
	  /*Update position of proteins in c2*/
	  translate_int(p2, c2, ind_com, bases, drev);
	  
	}//Done performing (or stopping) association
		
	
	//  cout <<"ind_com crds: "<<ind_com[c1].xcom<<' '<<ind_com[c1].ycom<<' '<<ind_com[c1].zcom<<endl;
	delete[] v0;
	delete[] v2;
	delete[] u;
	delete[] M;
	delete[] dtrans;
	delete[] drev;

	return cancel;
}
