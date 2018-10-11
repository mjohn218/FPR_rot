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

int associate_ap2_clath_sigma(int p1, int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom){
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//NOTE: This function is was writing by Osman and Basilio. However, this is a modification of associate_freelegPBCELL.cpp fuction///
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int prod = Rlist[mu][2];
	int iind = ihome[i1];
	int iind2 = ihome[i2];
	double Nv1;

	//get values of original associating complexes
	int c1 = bases[p1].mycomplex;
	int s1 = ind_com[c1].mysize;
	int c2 = bases[p2].mycomplex;
	int s2 = ind_com[c2].mysize;
	int newsize = s1 + s2;
	double trot;

	double *v1 = new double[3];
	double *v2 = new double[3];
	double *M = new double[9];
	double *Mneg = new double[9]; //same axis, negative angle
	double *u = new double[3];
	double *dtrans = new double[3];
	double *drev = new double[3];
	double *dtrans2 = new double[3];
	double *drev2 = new double[3];
	double cmpx=0.0, cmpy=0.0, cmpz=0.0;

	/*Now move protein*/
	double Dxtot = ind_com[c1].Dx + ind_com[c2].Dx;
	double Dytot = ind_com[c1].Dy + ind_com[c2].Dy;

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


	double theta;

	///// THIS WAS WORKING CODE COULD BE USED
//	//distance to move to place interfaces, along vector v
	dtrans[0] = -dx * ind_com[c1].Dx / Dxtot;//*(bindrad / R1 - 1.0);//place at sigma!
	dtrans[1] = -dy * ind_com[c1].Dy / Dytot;//*(bindrad / R1 - 1.0);//place at sigma!
	dtrans[2] = -dz * ind_com[c1].Dz / Dztot;//*(bindrad / R1 - 1.0);//place at sigma!

	drev[0] = +dx * ind_com[c2].Dx / Dxtot;//*(bindrad / R1 - 1.0);//place at sigma!
	drev[1] = +dy * ind_com[c2].Dy / Dytot;//*(bindrad / R1 - 1.0);//place at sigma!
	drev[2] = +dz * ind_com[c2].Dz / Dztot;//*(bindrad / R1 - 1.0);//place at sigma!
	if(bases[p1].protype == plist.pclath){//??? WHY JUST P1==? WHY NOT P2? It seems that if + ifelse (below) statement ensure that P1 or P2 are clath, if clath, change the center of mass????
	  
	  /*Uses the -3 to access the planar clathrin leg, rather than the out of plane clathrin-ap2 binding site*/
	  /////////////////////////
	  
	  cmpx = bases[p1].x[iind-3]-bases[p1].xcom;
	  cmpy = bases[p1].y[iind-3]-bases[p1].ycom;
	  cmpz = bases[p1].z[iind-3]-bases[p1].zcom;
	  //cout <<"In associate ap2_clath: iind of clathrin (p1): "<<iind<<" cmpx, y, z: "<<cmpx<<' '<<cmpy<<' '<<cmpz<<endl;
	  //BASILIO: Which values take iind? I Understand that iind is interfaces...so 1,2,3,4,5,6? what happen if iind take 1?
	  v1[0] = bases[p1].xcom+cmpx/2 - bases[p1].x[iind];
	  v1[1] = bases[p1].ycom+cmpy/2 - bases[p1].y[iind];
	  v1[2] = bases[p1].zcom+cmpz/2 - bases[p1].z[iind];
	  
	  
	  v2[0] = bases[p2].xcom - bases[p2].x[iind2];
	  v2[1] = bases[p2].ycom - bases[p2].y[iind2];
	  v2[2] = bases[p2].zcom - bases[p2].z[iind2];
	  
	  
	}else{
	  
	  
	  cmpx = bases[p2].x[iind2-3]-bases[p2].xcom;
	  cmpy = bases[p2].y[iind2-3]-bases[p2].ycom;
	  cmpz = bases[p2].z[iind2-3]-bases[p2].zcom;
	  
	  //cout <<"In associate ap2_clath: iind of clathrin (p2): "<<iind2<<" cmpx, y, z: "<<cmpx<<' '<<cmpy<<' '<<cmpz<<endl;
	  v2[0] = bases[p2].xcom+cmpx/2 - bases[p2].x[iind2];
	  v2[1] = bases[p2].ycom+cmpy/2 - bases[p2].y[iind2];
	  v2[2] = bases[p2].zcom+cmpz/2 - bases[p2].z[iind2];
	  
	  
	  v1[0] = bases[p1].xcom - bases[p1].x[iind];
	  v1[1] = bases[p1].ycom - bases[p1].y[iind];
	  v1[2] = bases[p1].zcom - bases[p1].z[iind];
	  
	}

	/*Calculate rotation matrix*/

	dotproduct(v2,v1, theta);
	//now calculate the axis of rotation
	crossproduct(v1, v2, u); //u is the UNIT vector of the rotation axis
	//double trot = 0.5 * (M_PI - theta); //to push them 180 apart           ///////////BASILIO MODIFIED THIS
	/*for rotating v1 around u the negative trot */
//	
	double Drztot=ind_com[c1].Drz+ind_com[c2].Drz;
	double trotneg, trotpos;
	if(Drztot<tol){
	  trotneg=(M_PI-theta)*0.5;
	  trotpos=trotneg;
	}else{
	  trotneg=(M_PI-theta)*ind_com[c1].Drz/Drztot;
	  trotpos=(M_PI-theta)*ind_com[c2].Drz/Drztot;
	}
	
	if(ind_com[c1].Dz!=0 && ind_com[c2].Dz!=0){	/*Update position of proteins in c2*/ // BOTH P1 AND P2 IN SOLUTION
	  
	  //trot = 0.5*(M_PI - theta); //to push ONE 180. In this case clath
	  
	  calc_Rmatrix(u, +trotpos,M);//c2
	  calc_Rmatrix(u, -trotneg,Mneg);//c1
	  rotate_and_translate_int(p1, c1, ind_com, bases, Mneg, dtrans, iind);
	  rotate_and_translate_int(p2, c2, ind_com, bases, M, drev, iind2); //BASILIO: I AM MODIFYING THIS FUNCION
	  
	  //		HERE WE ARE MISSING SIGMA DISPLACEMENT ALONG VECTOR V, TEST THIS
	  //ayudame(p1, p2, bases);
	  
	  /*The sign ensures the COM's get rotated
	    to the outside of interfaces final spot*/
	  
	  if(bases[p1].protype == plist.pclath){
	    
	    cmpx = bases[p1].x[iind-3]-bases[p1].xcom;
	    cmpy = bases[p1].y[iind-3]-bases[p1].ycom;
	    cmpz = bases[p1].z[iind-3]-bases[p1].zcom;
	    
	    
	    v1[0] = bases[p1].xcom+cmpx/2 - bases[p1].x[iind];
	    v1[1] = bases[p1].ycom+cmpy/2 - bases[p1].y[iind];
	    v1[2] = bases[p1].zcom+cmpz/2 - bases[p1].z[iind];
	    
	    
	  }else{
	    
	    v1[0] = bases[p1].xcom - bases[p1].x[iind];
	    v1[1] = bases[p1].ycom - bases[p1].y[iind];
	    v1[2] = bases[p1].zcom - bases[p1].z[iind];
	    
	    
	  }
	  Nv1 = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
	  //distance to move to place interfaces, along vector v
	  dtrans2[0] = v1[0] * ind_com[c1].Dx / Dxtot * bindrad/Nv1;
	  dtrans2[1] = v1[1] * ind_com[c1].Dy / Dytot * bindrad/Nv1;
	  dtrans2[2] = v1[2] * ind_com[c1].Dz / Dztot * bindrad/Nv1;
	  ; //should be zero
	  
	  drev2[0] = -v1[0] * ind_com[c2].Dx / Dxtot * bindrad/Nv1;
	  drev2[1] = -v1[1] * ind_com[c2].Dy / Dytot * bindrad/Nv1;
	  drev2[2] = -v1[2] * ind_com[c2].Dz / Dztot * bindrad/Nv1; //should be zero
	  
	  /*NOW WE ARE MOVING THE TWO PROTEINS TO A SEPARATION OF SIGMA*/
	  //just translate if larger complexes are coming together
	  translate_int(p1, c1, ind_com, bases, dtrans2);
	  
	  /*Same for complex2*********************************/
	  //just translate
	  translate_int(p2, c2, ind_com, bases, drev2);
	  //		ayudame(p1, p2, bases);
	  
	  
	  ///COMENTED BY BASILIO 11-07-2016
	}else if(ind_com[c1].Dz==0 && ind_com[c2].Dz==0){ // everybody's on the membrane ///???? I BELIEVE IT HAVE TO BE bases[p1].Dz==0, NOT !=0
	  /*No rotation or displacement in z*/
	  dtrans[2]=0;
	  drev[2]=0;
	  
	  calc_Rmatrix(u, 0, M);//c2
	  calc_Rmatrix(u, 0,Mneg);//c1
	  rotate_and_translate_int(p1, c1, ind_com, bases, Mneg, dtrans, iind);
	  rotate_and_translate_int(p2, c2, ind_com, bases, M, drev, iind2); //BASILIO: I AM MODIFYING THIS FUNCION
	  
	  // ???? I BELIEVE THAT NEXT SECTION SHOULD BE MODIFIED TO FIX PROBLEM CLATH + PIP2-AP2
	}else{//Clathrin in solution
	  ////what about clat-ap2 in solution interact with PI2?
	  //ADDED BY BASILIO and osman
	  //	  double trot = (M_PI - theta); //0.5 * (M_PI - theta);  //to push them 180 apart
	  /*for rotating v1 around u the negative trot */
	  //???MODIFY THIS JUST FOR CLAT SOLVE THE PROBLEM. I mean the angle MPI-THETA JUST FOR CLATH
	  
	  calc_Rmatrix(u, +trotpos,M);//c2
	  calc_Rmatrix(u, -trotneg,Mneg);//c1
	  
	  /*Update the position of proteins in complex one*/
	  // if(ind_com[c1].Dz!=0){//p1 is an clath molecule in the solution only rotate clathrin
	  rotate_and_translate_int(p1, c1, ind_com, bases, Mneg, dtrans, iind);
	  // }else{//if p2 is a clathrin molecule in the solution only rotate clathrin
	  /*Update position of proteins in c2*/
	  rotate_and_translate_int(p2, c2, ind_com, bases, M, drev, iind2); //BASILIO: I AM MODIFYING THIS FUNCION
	  //}
	  
	  ///////////////////////////////////////////////////////////////////////////
	  dx = -bases[p2].x[iind2] + bases[p1].x[iind];
	  dy = -bases[p2].y[iind2] + bases[p1].y[iind];
	  dz = -bases[p2].z[iind2] + bases[p1].z[iind];
	  //BAISLIO: DELTAX
	  
	  //boundary condition, inside box
	  
	  //distance to move to place interfaces, along vector v
	  dtrans2[0] = -dx * ind_com[c1].Dx / Dxtot;
	  dtrans2[1] = -dy * ind_com[c1].Dy / Dytot;
	  dtrans2[2] = -dz * ind_com[c1].Dz / Dztot;
	  
	  drev2[0] = +dx * ind_com[c2].Dx / Dxtot;
	  drev2[1] = +dy * ind_com[c2].Dy / Dytot;
	  drev2[2] = +dz * ind_com[c2].Dz / Dztot;
	  
	  translate_int(p1, c1, ind_com, bases, dtrans2);
	  translate_int(p2, c2, ind_com, bases, drev2);
	  dtrans2[0]=0;
	  dtrans2[1]=0;
	  dtrans2[2]=bindrad;
	  
	  if(ind_com[c1].Dz!=0){///if p1 isnot Ap2(IS CLATHRIN) do this
	    translate_int(p1, c1, ind_com, bases, dtrans2);
	  }else{//if p2 isnot Ap2 do this
		/*Update position of proteins in c2*/
	    translate_int(p2, c2, ind_com, bases, dtrans2);
	  }
	  ///////////////////////////////////////////////////////////////////////////
	  //		ayudame(p1, p2, bases);
	}
		
	int cancel;
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
	  /*reverse initial rotation and translation*/
	  for(int i=0;i<3;i++){
	    dtrans[i]*=-1;
	    drev[i]*=-1;
	  }
	  calc_Rmatrix(u, -trotpos,M);//c2
	  calc_Rmatrix(u, +trotneg,Mneg);//c1
	  rotate_and_translate_int(p1, c1, ind_com, bases, Mneg, dtrans, iind);
	  rotate_and_translate_int(p2, c2, ind_com, bases, M, drev, iind2); //BASILIO: I AM MODIFYING THIS FUNCION
	  
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
	delete[] dtrans2;
	delete[] drev2;
	return cancel;
}
