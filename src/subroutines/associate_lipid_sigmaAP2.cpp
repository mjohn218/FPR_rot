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

int associate_lipid_sigmaAP2(int p1, int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom){

  /*THIS ROUTINE IS ONLY CALLED WHEN BINDING TO A LIPID. GIVE LIPID ORITENTATION A VECTOR [0,0,1], SO IT DOESN'T NEED A SEPARATE INTERFACE.*/
  /*March 2017, put them at sigma when they associate!*/
	int prod = Rlist[mu][2];
	int iind = ihome[i1];
	int iind2 = ihome[i2];

	//get values of original associating complexes
	int c1 = bases[p1].mycomplex;
	int c2 = bases[p2].mycomplex;

	int cancel;

	double *v1 = new double[3];
	double *v2 = new double[3];
	double *M = new double[9];
	double *Mneg = new double[9]; //same axis, negative angle
	double *u = new double[3];
	double *dtrans = new double[3];
	double *drev = new double[3];
	double *dtrans1 = new double[3];
	double *drev1 = new double[3];

	/*Now move protein*/
	double Dxtot = ind_com[c1].Dx + ind_com[c2].Dx;
	double Dytot = ind_com[c1].Dy + ind_com[c2].Dy;
	//Dxtot:total diffusion constant in x

	double Dztot = ind_com[c1].Dz + ind_com[c2].Dz;
	double tol = 1E-16;
	if (Dztot < tol)
		Dztot = 1; //otherwise you divide by zero

	//distance between the associating proteins interfaces!
	/*The sign needs to be switched to ensure the COM's get rotated
	 to the outside of interfaces final spot*/
	double dx = -bases[p2].x[iind2] + bases[p1].x[iind];
	double dy = -bases[p2].y[iind2] + bases[p1].y[iind];
	double dz = -bases[p2].z[iind2] + bases[p1].z[iind];
	//BAISLIO: DELTAX
	
	/*Put at contact!*/
	double R2 = dx*dx+dy*dy+dz*dz;
	double R1 = sqrt(R2);
	

	//boundary condition, inside box

	double theta;

	//distance to move to place interfaces, along vector v//
	dtrans1[0] = -dx * ind_com[c1].Dx / Dxtot;
	dtrans1[1] = -dy * ind_com[c1].Dy / Dytot;
	dtrans1[2] = -dz * ind_com[c1].Dz / Dztot;
	
	drev1[0] = +dx * ind_com[c2].Dx / Dxtot;
	drev1[1] = +dy * ind_com[c2].Dy / Dytot;
	drev1[2] = +dz * ind_com[c2].Dz / Dztot;


	if(bases[p1].Dz==0){
	  //BASILIO: VECTOR V1 COORDINATES
	  v1[0] = 0;
	  v1[1] = 0;
	  v1[2] = -1;
	  
	  v2[0] = bases[p2].xcom - bases[p2].x[iind2];
	  v2[1] = bases[p2].ycom - bases[p2].y[iind2];
	  v2[2] = bases[p2].zcom - bases[p2].z[iind2];
	  
	  //BASILIO:TIL HERE
	}else{
	  //BASILIO:VECOR V2 COORDINATES
	  v1[0] = bases[p1].xcom - bases[p1].x[iind];
	  v1[1] = bases[p1].ycom - bases[p1].y[iind];
	  v1[2] = bases[p1].zcom - bases[p1].z[iind];
	  
	  //BASILIO:TIL HERE
	  v2[0] = 0;
	  v2[1] = 0;
	  v2[2] = -1;
	  
	}

	/*Calculate rotation matrix*/

	dotproduct(v2,v1, theta);
	double Drztot=ind_com[c1].Drz+ind_com[c2].Drz;
	double trotneg, trotpos;
	if(Drztot<tol){
	  trotneg=(M_PI-theta)*0.5;
	  trotpos=trotneg;
	  
	}else{

	  trotneg=(M_PI-theta)*ind_com[c1].Drz/Drztot;
	  trotpos=(M_PI-theta)*ind_com[c2].Drz/Drztot;
	}
	
	//ayudame(p1, p2, bases);
	//now calculate the axis of rotation
	crossproduct(v1, v2, u); //u is the UNIT vector of the rotation axis
	//double trot = (M_PI - theta); //0.5 * (M_PI - theta);  //to push them 180 apart
	/*for rotating v1 around u the negative trot */
	//???MODIFY THIS JUST FOR CLAT SOLVE THE PROBLEM. I mean the angle MPI-THETA JUST FOR CLATH

	calc_Rmatrix(u, +trotpos,M);//c2
	calc_Rmatrix(u, -trotneg, Mneg);//c1
	
	/*Update the position of proteins in complex one*/
	//if(bases[p1].zcom!=-(plist.zboxl/2)){///if p1 isnot pip2 do this
	rotate_and_translate_int(p1, c1, ind_com, bases, Mneg, dtrans1, iind);
	//}else{//if p2 isnot pip2 do this
	/*Update position of proteins in c2*/
	rotate_and_translate_int(p2, c2, ind_com, bases, M, drev1, iind2); //BASILIO: I AM MODIFYING THIS FUNCION
	//}
	///???Do we need the above section? (I mean the above 4lines?)

///////////////////////////////////////////////////////////////////////////
	/*Needed this below to put them at sigma in z.*/
// dx = bases[p2].xcom - bases[p1].xcom;
	// dy = bases[p2].ycom - bases[p1].ycom;
	// dz = bases[p2].zcom - bases[p1].zcom;
	// R2 = dx*dx+dy*dy+dz*dz;
	// R1 = sqrt(R2);

	cout <<"LIPID to AP2, after rot and translate, current sep: "<<dx<<' '<<dy<<' '<<dz<<endl;
	//BAISLIO: DELTAX
	double phi=rand_gsl()*2.0*M_PI;//between 0 and 2pi
	dx=bindrad*cos(phi);
	dy=bindrad*sin(phi);
	dz=0;
	/*Move APART by bindrad, put them adjacent, i.e., move in x,y, not in z!*/	
	dtrans[0] = -dx * ind_com[c1].Dx / Dxtot;//*bindrad / R1;//place at sigma!
	dtrans[1] = -dy * ind_com[c1].Dy / Dytot;//*bindrad / R1;//place at sigma!
	dtrans[2] = -dz * ind_com[c1].Dz / Dztot;//*bindrad / R1;//place at sigma!

	drev[0]=+dx * ind_com[c2].Dx / Dxtot;//*bindrad / R1;//place at sigma!
	drev[1]=+dy * ind_com[c2].Dy / Dytot;//*bindrad / R1;//place at sigma!
	drev[2]=+dz * ind_com[c2].Dz / Dztot;//*bindrad / R1;//place at sigma!
	
	translate_int(p1, c1, ind_com, bases, dtrans);
	translate_int(p2, c2, ind_com, bases, drev);
	
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
	  cout << "Unbind, clathrins crashed together ! " << endl;
	  /*Reverse translation and rotation*/
	  for(int i=0;i<3;i++){
	    dtrans[i]*=-1;
	    drev[i]*=-1;
	    dtrans1[i]*=-1;
	    drev1[i]*=-1;
	  }
	  /*reverse last translation*/
	  translate_int(p1, c1, ind_com, bases, dtrans);
	  translate_int(p2, c2, ind_com, bases, drev);
	  /*Reverse rotation*/
	  calc_Rmatrix(u, -trotpos,M);//c2
	  calc_Rmatrix(u, +trotneg, Mneg);//c1
	  
	  rotate_and_translate_int(p1, c1, ind_com, bases, Mneg, dtrans1, iind);
	  rotate_and_translate_int(p2, c2, ind_com, bases, M, drev1, iind2);
	  /*Update complex com's, just in case, and check it is inside the box.*/
	  update_one_com_only(c1, ind_com, bases);
	  update_one_com_only(c2, ind_com, bases);
	  reflect_complex_rad_rot(p1, bases, ind_com, plist);
	  reflect_complex_rad_rot(p2, bases, ind_com, plist);
	  
	
	}//Done performing (or stopping) association
	
	cout << "final complex com: " << ind_com[c1].xcom << ' ' << ind_com[c1].ycom << ' ' << ind_com[c1].zcom << " radius: " << ind_com[c1].radR << " Dr: " << ind_com[c1].Drx << " Dtrans: " << ind_com[c1].Dx << endl;

//	delete[] u1;
//	delete[] u2;
	delete[] v1;
	delete[] v2;
	delete[] u;
	delete[] M;
	delete[] Mneg;
	delete[] dtrans;
	delete[] drev;
	return cancel;
}
