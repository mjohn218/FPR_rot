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

int associate_2dclat_3dclat_lipid_sigma(int p1, int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//NOTE: This function is was writing by Basilio. However, this is a modification of associate_freelegPBCELL.cpp fuction///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
	double *uu=new double[3];
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
	
	double theta;

	//distance to move to place interfaces, along vector v
	dtrans1[0] = -dx * ind_com[c1].Dx / Dxtot;//*(bindrad / R1 - 1.0);//place at sigma!
	dtrans1[1] = -dy * ind_com[c1].Dy / Dytot;//*(bindrad / R1 - 1.0);//place at sigma!
	dtrans1[2] = -dz * ind_com[c1].Dz / Dztot;//*(bindrad / R1 - 1.0);//place at sigma!

	drev1[0] = +dx * ind_com[c2].Dx / Dxtot;//*(bindrad / R1 - 1.0);//place at sigma!
	drev1[1] = +dy * ind_com[c2].Dy / Dytot;//*(bindrad / R1 - 1.0);//place at sigma!
	drev1[2] = +dz * ind_com[c2].Dz / Dztot;


	//BASILIO: VECTOR V1 COORDINATES
	v1[0] = bases[p1].xcom - bases[p1].x[iind];
	v1[1] = bases[p1].ycom - bases[p1].y[iind];
	v1[2] = bases[p1].zcom - bases[p1].z[iind];
	//BASILIO:TIL HERE

	//BASILIO:VECOR V2 COORDINATES
	v2[0] = bases[p2].xcom - bases[p2].x[iind2];
	v2[1] = bases[p2].ycom - bases[p2].y[iind2];
	v2[2] = bases[p2].zcom - bases[p2].z[iind2];
	//BASILIO:TIL HERE
	/*Calculate rotation matrix*/

	dotproduct(v2,v1, theta);
	//now calculate the axis of rotation
	crossproduct(v1, v2, u); //u is the UNIT vector of the rotation axis
	double Drztot=ind_com[c1].Drz+ind_com[c2].Drz;
	double trotneg, trotpos;
	if(Drztot<tol){
	  trotneg=(M_PI-theta)*0.5;
	  trotpos=trotneg;
	}else{
	  trotneg=(M_PI-theta)*ind_com[c1].Drz/Drztot;
	  trotpos=(M_PI-theta)*ind_com[c2].Drz/Drztot;
	}
	
	//double trot = (M_PI - theta); //0.5 * (M_PI - theta);  //to push them 180 apart
	/*for rotating v1 around u the negative trot */
	//???MODIFY THIS JUST FOR CLAT SOLVE THE PROBLEM. I mean the angle MPI-THETA JUST FOR CLATH
	
	calc_Rmatrix(u, +trotpos,M); //for c2
	calc_Rmatrix(u, -trotneg, Mneg);//for c1
	
	/*WITH DRZ WEIGHTING WE SHOULD NOT NEED THE BELOW IF STATEMENTS*/
	/*Update the position of proteins in complex one*/
	if(ind_com[bases[p1].mycomplex].Dz!=0){
	  rotate_and_translate_int(p1, c1, ind_com, bases, Mneg, dtrans1, iind);
	  translate_int(p2, c2, ind_com, bases, drev1);//move to contact
	}else{
	  rotate_and_translate_int(p2, c2, ind_com, bases, M, drev1, iind2);
	  translate_int(p1, c1, ind_com, bases, dtrans1);//move to contact
	}
	
	///???Do we need the above section? (I mean the above 4lines?)
	/*nOW PUSH THE INTERFACES TO SEPARATION SIGMA, ALONG THE COM-COM VECTOR*/
	dx = bases[p2].xcom - bases[p1].xcom;
	dy = bases[p2].ycom - bases[p1].ycom;
	dz = bases[p2].zcom - bases[p1].zcom;
	R2 = dx*dx+dy*dy+dz*dz;
	R1 = sqrt(R2);
	
	//We want to move them apart from each other, along this vector.
	dtrans[0] = -dx * ind_com[c1].Dx / Dxtot*bindrad / R1;//place at sigma!
	dtrans[1] = -dy * ind_com[c1].Dy / Dytot*bindrad / R1;//place at sigma!
	dtrans[2] = -dz * ind_com[c1].Dz / Dztot*bindrad / R1;//place at sigma!

	drev[0] = +dx * ind_com[c2].Dx / Dxtot*bindrad / R1 ;//place at sigma!
	drev[1] = +dy * ind_com[c2].Dy / Dytot*bindrad / R1 ;//place at sigma!
	drev[2] = +dz * ind_com[c2].Dz / Dztot*bindrad/R1;

	
	translate_int(p1, c1, ind_com, bases, dtrans);//move to contact
	translate_int(p2, c2, ind_com, bases, drev);
	
	/*Rotate to proper orientation*/
	int ind3;
	double *u1 = new double[3];
	double *u2 = new double[3];
	  ////////////////////////////////////////////////////////////////////////////////////
	  ///////////////////////////// MODIFY BY BASILIO ////////////////////////////////////
	// This section use to use iind values of 1,2 and 3. However, now we are using 0,1,2. Previous version below
	if (iind == 0)
		ind3 = 2;
	else
		ind3 = iind - 1;
	/////////////////////////////
//	PREVIOUS VERSION

//	if (iind == 1)
//		ind3 = 3;
//	else
//		ind3 = iind - 1;
/////////////////////////////////////////////
	  ///////////////////// MODIFICATION END HERE /////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////
	v2[0] = bases[p1].x[iind] - bases[p1].xcom;
	v2[1] = bases[p1].y[iind] - bases[p1].ycom;
	v2[2] = bases[p1].z[iind] - bases[p1].zcom;


	v1[0] = bases[p1].x[ind3] - bases[p1].xcom;
	v1[1] = bases[p1].y[ind3] - bases[p1].ycom;
	v1[2] = bases[p1].z[ind3] - bases[p1].zcom;


	/*2 cross 1 will be normal up*/
	crossproduct(v2, v1, u1);
	cout << "normal 1: " << u1[0] << ' ' << u1[1] << ' ' << u1[2] << endl;
	v2[0] = bases[p2].x[2] - bases[p2].xcom;
	v2[1] = bases[p2].y[2] - bases[p2].ycom;
	v2[2] = bases[p2].z[2] - bases[p2].zcom;


	v1[0] = bases[p2].x[1] - bases[p2].xcom;
	v1[1] = bases[p2].y[1] - bases[p2].ycom;
	v1[2] = bases[p2].z[1] - bases[p2].zcom;


	crossproduct(v2, v1, u2);
	cout << "normal 2: " << u2[0] << ' ' << u2[1] << ' ' << u2[2] << endl;
	dotproduct(u2,u1, theta);
	double trotneg1, trotpos1;
	if(Drztot<tol){
	  trotneg1=(theta)*0.5;
	  trotpos1=trotneg1;
	}else{
	  trotneg1=(theta)*ind_com[c2].Drz/Drztot;
	  trotpos1=(theta)*ind_com[c1].Drz/Drztot;
	}

	//	trot = theta;
	double *M2 = new double[9];
	double *M2neg = new double[9];
	/*p1 rotate by positive theta/2, p2 rotate by negative theta/2
	 rotate around the bound leg segment
	 */
	crossproduct(u1, u2, uu);
	calc_Rmatrix(uu, +trotpos1, M2);//c1
	calc_Rmatrix(uu, -trotneg1, M2neg);//c2
	
	if(ind_com[bases[p1].mycomplex].Dz!=0 && ind_com[bases[p2].mycomplex].Dz==0){//p1 is an clath molecule in the solution only rotate clathrin
	  rotate_only(p1, c1, ind_com, bases, M2);
	}else{//if p2 is a clathrin molecule in the solution only rotate clathrin
	/*Update position of proteins in c2*/
	  rotate_only(p2, c2, ind_com, bases, M2neg);
	}
	
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
	  
	  /*CANCEL=1: IN THIS CASE CLATHRINS CRASHED TOGETHER
	    IN THIS CASE< CLATHRINS CRASHED TOGETHER
	    un-bind them, allow to diffuse apart*/
	  cout << "Unbind clathrins, crashed together ! " << endl;
	  /*reverse previous rotation*/
	  calc_Rmatrix(uu, -trotpos1, M2);//c1
	  calc_Rmatrix(uu, +trotneg1, M2neg);//c2
	  
	  if(ind_com[bases[p1].mycomplex].Dz!=0 && ind_com[bases[p2].mycomplex].Dz==0){//p1 is an clath molecule in the solution only rotate clathrin
	    rotate_only(p1, c1, ind_com, bases, M2);
	  }else{//if p2 is a clathrin molecule in the solution only rotate clathrin
	    /*Update position of proteins in c2*/
	    rotate_only(p2, c2, ind_com, bases, M2neg);
	  }
	  for(int i=0;i<3;i++){
	    dtrans[i]*=-1;
	    drev[i]*=-1;
	    dtrans1[i]*=-1;
	    drev1[i]*=-1;
	  }
	  /*reverse translation*/	  
	  translate_int(p1, c1, ind_com, bases, dtrans);//move to contact
	  translate_int(p2, c2, ind_com, bases, drev);
	  
	  /*Reverse initial rotation*/
	  calc_Rmatrix(u, -trotpos,M); //for c2
	  calc_Rmatrix(u, +trotneg, Mneg);//for c1
	  
	  /*WITH DRZ WEIGHTING WE SHOULD NOT NEED THE BELOW IF STATEMENTS*/
	  /*Update the position of proteins in complex one*/
	  if(ind_com[bases[p1].mycomplex].Dz!=0){
	    rotate_and_translate_int(p1, c1, ind_com, bases, Mneg, dtrans1, iind);
	    translate_int(p2, c2, ind_com, bases, drev1);//move to contact
	  }else{
	    rotate_and_translate_int(p2, c2, ind_com, bases, M, drev1, iind2);
	    translate_int(p1, c1, ind_com, bases, dtrans1);//move to contact
	  }
	  /*Update complex com's, just in case, and check it is inside the box.*/
	  update_one_com_only(c1, ind_com, bases);
	  update_one_com_only(c2, ind_com, bases);
	  reflect_complex_rad_rot(p1, bases, ind_com, plist);
	  reflect_complex_rad_rot(p2, bases, ind_com, plist);
	  
	}
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
