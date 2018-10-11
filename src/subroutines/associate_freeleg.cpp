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

int associate_freeleg(int p1, int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom)
{

	int prod = Rlist[mu][2];
	int iind = ihome[i1];
	int iind2 = ihome[i2];


	//get values of original associating complexes
	int c1 = bases[p1].mycomplex;
	int c2 = bases[p2].mycomplex;

	double *v1 = new double[3];
	double *v2 = new double[3];
	double *M = new double[9];
	double *Mneg = new double[9]; //same axis, negative angle
	double *u = new double[3];
	double *uu = new double[3];
	double *dtrans = new double[3];
	double *drev = new double[3];
	int cancel;


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
//BASILIO: This(Above) is working for all the intefaces? So...iind is all the interfaces? because I do not see a loop?

	/*Put at contact!*/
	double R2 = dx*dx+dy*dy+dz*dz;
	double R1 = sqrt(R2);
	
	double theta;

	//distance to move to place interfaces, along vector v
	dtrans[0] = -dx * ind_com[c1].Dx / Dxtot;
	dtrans[1] = -dy * ind_com[c1].Dy / Dytot;
	dtrans[2] = -dz * ind_com[c1].Dz / Dztot;

	drev[0] = +dx * ind_com[c2].Dx / Dxtot;
	drev[1] = +dy * ind_com[c2].Dy / Dytot;
	drev[2] = +dz * ind_com[c2].Dz / Dztot;

	v1[0] = bases[p1].xcom - bases[p1].x[iind];
	v1[1] = bases[p1].ycom - bases[p1].y[iind];
	v1[2] = bases[p1].zcom - bases[p1].z[iind];

	v2[0] = bases[p2].xcom - bases[p2].x[iind2];
	v2[1] = bases[p2].ycom - bases[p2].y[iind2];
	v2[2] = bases[p2].zcom - bases[p2].z[iind2];

	/*Calculate rotation matrix*/

	/*If v1 and v2 are antiparallel, no need to rotate (theta=pi), if they are parallel (theta=0)
	  they need to rotate by pi, but cross product will not give unique vector for rotation axis.
	 */
	dotproduct(v2,v1, theta);
	
	if(M_PI-theta<1E-8 || theta<1E-8){
	  //no rotation, already aligned for pi, parallel for 0. No cross prod in either case. Use z-axis
	  uu[0]=0;
	  uu[1]=0;
	  uu[2]=1;
	}else{
	  //now calculate the axis of rotation
	  crossproduct(v1, v2, uu); //u is the UNIT vector of the rotation axis
	}
	
	if(ind_com[c1].Dz<tol && ind_com[c2].Dz<tol){
	  /*Make sure that uu is a normal facing only in z, so no rotation occurs around x or y!*/
	  if(abs(uu[2])!=1){
	    cout.precision(12);
	    cout <<"WARNING: On membrane, normal vector: "<<uu[0]<<' '<<uu[1]<<' '<<uu[2]<<endl;
	    cout <<" But it should point only in z "<<endl;
	    cout <<" v1: "<<v1[0]<<' '<<v1[1]<<' '<<v1[2]<<endl;
	    cout <<" v2: "<<v2[0]<<' '<<v2[1]<<' '<<v2[2]<<endl;
	    cout <<"theta:" <<theta<<endl;
	    cout <<"Correcting..."<<endl;
	    uu[0]=0;
	    uu[1]=0;
	    uu[2]=round(uu[2]);//could be + or -1

	  }
	}

	double Drztot=ind_com[c1].Drz+ind_com[c2].Drz;
	double trotneg, trotpos;
	double trotneg1, trotpos1;
	if(Drztot<tol){
	  trotneg1=(M_PI-theta)*0.5;
	  trotpos1=trotneg1;
	  
	}else{

	  trotneg1=(M_PI-theta)*ind_com[c1].Drz/Drztot;
	  trotpos1=(M_PI-theta)*ind_com[c2].Drz/Drztot;
	}
	//double trot = 0.5 * (M_PI - theta); //to push them 180 apart
	/*for rotating v1 around u the negative trot */

	calc_Rmatrix(uu, +trotpos1,M);//for c2
	calc_Rmatrix(uu, -trotneg1, Mneg);//for c1
	
	/*Update the position of proteins in complex one*/
	rotate_and_translate_int(p1, c1, ind_com, bases, Mneg, dtrans, iind);

	/*Update position of proteins in c2*/
	rotate_and_translate_int(p2, c2, ind_com, bases, M, drev, iind2);

	
	/*NOW PUSH THE INTERFACES TO SEPARATION SIGMA, ALONG THE COM-COM VECTOR*/
	dx = bases[p2].xcom - bases[p1].xcom;
	dy = bases[p2].ycom - bases[p1].ycom;
	dz = bases[p2].zcom - bases[p1].zcom;
	R2 = dx*dx+dy*dy+dz*dz;
	R1 = sqrt(R2);
	//We want to move them APART from each other, along this vector.
	double *dtrans2=new double[3];
	double *drev2=new double[3];
	dtrans2[0] = -dx * ind_com[c1].Dx / Dxtot*bindrad / R1;//place at sigma!
	dtrans2[1] = -dy * ind_com[c1].Dy / Dytot*bindrad / R1;//place at sigma!
	dtrans2[2] = -dz * ind_com[c1].Dz / Dztot*bindrad / R1;//place at sigma!

	//cout <<"to displace: "<<dtrans2[0]<<' '<<dtrans2[1]<<' '<<dtrans2[2]<<endl;
	drev2[0]=+dx * ind_com[c2].Dx / Dxtot*bindrad / R1;//place at sigma!
	drev2[1]=+dy * ind_com[c2].Dy / Dytot*bindrad / R1;//place at sigma!
	drev2[2]=+dz * ind_com[c2].Dz / Dztot*bindrad / R1;//place at sigma!
	
	translate_int(p1, c1, ind_com, bases, dtrans2);//move to contact
	translate_int(p2, c2, ind_com, bases, drev2);
	
	/*Rotate to proper orientation*///??? what this mean???

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
	///////////////////////////
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
	//cout << "normal 1: " << u1[0] << ' ' << u1[1] << ' ' << u1[2] << endl;
	if(ind_com[c1].Dz<tol && ind_com[c2].Dz<tol){
	  /*Make sure that uu is a normal facing only in z, so no rotation occurs around x or y!*/
	  if(abs(u1[2])!=1){
	    cout.precision(12);
	    cout <<"WARNING: On membrane, normal vector: "<<u1[0]<<' '<<u1[1]<<' '<<u1[2]<<endl;
	    cout <<" But it should point only in z "<<endl;
	    cout <<" v1: "<<v1[0]<<' '<<v1[1]<<' '<<v1[2]<<endl;
	    cout <<" v2: "<<v2[0]<<' '<<v2[1]<<' '<<v2[2]<<endl;
	    cout <<"theta:" <<theta<<endl;
	    cout <<"Correcting..."<<endl;
	    u1[0]=0;
	    u1[1]=0;
	    u1[2]=round(u1[2]);//could be + or - 1

	  }
	}

	v2[0] = bases[p2].x[2] - bases[p2].xcom;
	v2[1] = bases[p2].y[2] - bases[p2].ycom;
	v2[2] = bases[p2].z[2] - bases[p2].zcom;


	v1[0] = bases[p2].x[1] - bases[p2].xcom;
	v1[1] = bases[p2].y[1] - bases[p2].ycom;
	v1[2] = bases[p2].z[1] - bases[p2].zcom;


	crossproduct(v2, v1, u2);
	//cout << "normal 2: " << u2[0] << ' ' << u2[1] << ' ' << u2[2] << endl;

	if(ind_com[c1].Dz<tol && ind_com[c2].Dz<tol){
	  /*Make sure that uu is a normal facing only in z, so no rotation occurs around x or y!*/
	  if(abs(u2[2])!=1){
	    cout.precision(12);
	    cout <<"WARNING: On membrane, normal vector: "<<u2[0]<<' '<<u2[1]<<' '<<u2[2]<<endl;
	    cout <<" But it should point only in z "<<endl;
	    cout <<" v1: "<<v1[0]<<' '<<v1[1]<<' '<<v1[2]<<endl;
	    cout <<" v2: "<<v2[0]<<' '<<v2[1]<<' '<<v2[2]<<endl;
	    cout <<"theta:" <<theta<<endl;
	    cout <<"Correcting..."<<endl;
	    u2[0]=0;
	    u2[1]=0;
	    u2[2]=round(u2[2]);

	  }
	}
	
	dotproduct(u2,u1, theta);
	if(Drztot<tol){
	  trotneg=(theta)*0.5;
	  trotpos=trotneg;
	  
	}else{

	  trotneg=(theta)*ind_com[c2].Drz/Drztot;
	  trotpos=(theta)*ind_com[c1].Drz/Drztot;
	}
	

	double *M2 = new double[9];
	double *M2neg = new double[9];
	/*p1 rotate by positive theta/2, p2 rotate by negative theta/2
	 rotate around the bound leg segment
	 */
	crossproduct(u1, u2, u);
	calc_Rmatrix(u, +trotpos, M2);//c1
	calc_Rmatrix(u, -trotneg, M2neg);//c2
	//cout <<"Rotate c1 by: "<<trotpos<<" and c2 by : "<<-trotneg<<endl;
	
	rotate_only(p1, c1, ind_com, bases, M2);
	rotate_only(p2, c2, ind_com, bases, M2neg);
	
	/*Determine if you crashed two clathrins together
	  unless they are in the same complex, in which case you are closing a loop.
	*/
	//cout <<"Current positions of clathrins, before checking overlap: "<<endl;
	update_one_com_only(c1,ind_com, bases);
	update_one_com_only(c2,ind_com, bases);

	//write_crds_complex(c1, ind_com, bases);
	//write_crds_complex(c2, ind_com, bases);
	
	cancel=0;
	if (c1 != c2)
	  cancel = measure_overlap(c1, c2, ind_com, bases, plist.pclath);//Not closing a loop, so check for overlap clashes
	if(cancel==0)
	  cancel=check_box_span(p1, p2, bases, ind_com, plist);//if complex is too big to fit in the box, reverse the association!

	if(cancel==0)
	  update_bound_proteins( p1,  p2,  i1,  i2,  bases,  ind_com,  plist,  ncrosscom, prod,  iind, iind2);
	else {
	  /*IN THIS CASE< CLATHRINS CRASHED TOGETHER
	    un-bind them, allow to diffuse apart*/
	  if(cancel==-1)cout <<"Spanned box, unbind clathrins "<<endl;
	  else
	    cout << "Unbind clathrins, crashed together ! " << endl;
	  /*reverse previous rotation*/
	  calc_Rmatrix(u, -trotpos, M2);//c1
	  calc_Rmatrix(u, +trotneg, M2neg);//c2
	
	  rotate_only(p1, c1, ind_com, bases, M2);
	  rotate_only(p2, c2, ind_com, bases, M2neg);
	  
	  /*reverse translation*/
	  for(int i=0;i<3;i++){
	    dtrans[i]*=-1;
	    drev[i]*=-1;
	    dtrans2[i]*=-1;
	    drev2[i]*=-1;
	  }
	  translate_int(p1, c1, ind_com, bases, dtrans2);//move to contact
	  translate_int(p2, c2, ind_com, bases, drev2);
	  
	  calc_Rmatrix(uu, -trotpos1,M);//for c2
	  calc_Rmatrix(uu, +trotneg1, Mneg);//for c1
	
	  /*Reverse initial rotation and translation*/
	  rotate_and_translate_int(p1, c1, ind_com, bases, Mneg, dtrans, iind);
	  rotate_and_translate_int(p2, c2, ind_com, bases, M, drev, iind2);
	  /*Update complex com's, just in case, and check it is inside the box.*/
	  update_one_com_only(c1, ind_com, bases);
	  update_one_com_only(c2, ind_com, bases);
	  reflect_complex_rad_rot(p1, bases, ind_com, plist);
	  reflect_complex_rad_rot(p2, bases, ind_com, plist);
	  
	}

	cout << "final complex com: " << ind_com[c1].xcom << ' ' << ind_com[c1].ycom << ' ' << ind_com[c1].zcom << " radius: " << ind_com[c1].radR << " Dr: " << ind_com[c1].Drx << " Dtrans: " << ind_com[c1].Dx << endl;

	delete[] u1;
	delete[] uu;
	delete[] u2;
	delete[] v1;
	delete[] v2;
	delete[] u;
	delete[] M;
	delete[] Mneg;
	delete[] dtrans;
	delete[] drev;
	return cancel;
}
