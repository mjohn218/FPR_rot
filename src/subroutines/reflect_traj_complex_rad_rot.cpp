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

void reflect_traj_complex_rad_rot(int p1, Fullmol *bases, Complex *ind_com, Parms plist, double *traj, double *trajR, double *M) {
  /*This routine updated March 2017 to test if a large complex that spans the box could extend out in both directions
    if so, it attempts to correct for this by resampling the complex's translational and rotational updates. 
   */
  int i, j;
	int k = bases[p1].mycomplex;
	int s1 = ind_com[k].mysize;
	double xboxl=plist.xboxl;
	double yboxl=plist.yboxl;
	double zboxl=plist.zboxl;
	
	double xchg;
	double ychg;
	double zchg;
	double xtot = 0.0;
	double ytot = 0.0;
	double ztot = 0.0;
	int flagx = 0;
	int flagy = 0;
	double rad = ind_com[k].radR;
	int flagz = 0;
	double currx = ind_com[k].xcom + traj[0];
	double curry = ind_com[k].ycom + traj[1];
	double currz = ind_com[k].zcom + traj[2];
	double xtot0, ytot0, ztot0;
	double tol = 1E-11;	

	/*This is to test based on general size if it is close to boundaries, before doing detailed evaluation below.*/

	xchg = currx - xboxl / 2.0;
	if ((xchg + rad) > 0) 
	  flagx = 1;
	if ((xchg - rad) < -xboxl) 
	  flagx += 1;
	
	ychg = curry - yboxl / 2.0;
	if ((ychg + rad) > 0) 
	  flagy = 1;
	if ((ychg - rad) < -yboxl) 
	  flagy += 1;
	
	zchg = currz - zboxl / 2.0;
	if ((zchg + rad) > 0) 
	  flagz = 1;
	if ((zchg - rad) < -zboxl) 
		flagz += 1;
	

	int mp;

	/*Now evaluate all interfaces distance from boundaries.*/
	int recheck=0;
	double dx, dy, dz;
	double x0, y0, z0;
	double dzrot, dxrot, dyrot;
	double row[3];
	int flag;
	int fneg, fpos;
	double maxpos, maxneg;
	double poswall, negwall;
	if (flagx >0) {
	  flag = 0;
	  xtot = 0;
	  fneg=0;
	  fpos=0;
	  poswall=plist.xboxl;
	  negwall=plist.xboxl;
	  x0 = ind_com[k].xcom;
	  y0 = ind_com[k].ycom;
	  z0 = ind_com[k].zcom;
	  row[0] = M[0];
	  row[1] = M[1];
	  row[2] = M[2];
	  
	  /*these need to be what current positions
	    due to translation and rotation are*/
	  for (i = 0; i < s1; i++) {
	    mp = ind_com[k].plist[i];
	    
	    /*measure each protein COM to z plane*/
	    dx = bases[mp].xcom - x0;
	    dy = bases[mp].ycom - y0;
	    dz = bases[mp].zcom - z0;
	      
	    dxrot = row[0] * dx + row[1] * dy + row[2] * dz;
	    currx = x0 + traj[0] + dxrot;
	    if(xboxl/2.0-currx<poswall)
	      poswall=xboxl/2.0-currx;//shortest distance in x from the +wall
	    if(currx+xboxl/2.0<negwall)
	      negwall=currx+xboxl/2.0;//shortest distance in x from the -wall. >0 inside the box, <0 outside the box.
	    
	    xchg = currx - xboxl / 2.0;
	    if (xchg > 0) {
	      fpos=1;
	      flag = 1;
	      if (-xchg < xtot)
		xtot = -xchg;
	      
	    } else if (xchg < -xboxl) {
	      fneg=1;
	      flag = 1;
	      if (-(xchg + xboxl) > xtot)
		xtot = -(xchg + xboxl);
	      
	    }
	    /*measure each interface to z plane*/
	    for (j = 0; j < bases[mp].ninterface; j++) {
	      dx = bases[mp].x[j] - x0;
	      dy = bases[mp].y[j] - y0;
	      dz = bases[mp].z[j] - z0;
	      
	      dxrot = row[0] * dx + row[1] * dy + row[2] * dz;
	      currx = x0 + traj[0] + dxrot;
	      if(xboxl/2.0-currx<poswall)
		poswall=xboxl/2.0-currx;//shortest distance in x from the +wall
	      if(currx+xboxl/2.0<negwall)
		negwall=currx+xboxl/2.0;//shortest distance in x from the -wall. >0 inside the box, <0 outside the box.
	      
	      xchg = currx - xboxl / 2.0;
	      if (xchg > 0) {
		fpos=1;
		flag = 1;
		if (-xchg < xtot)
		  xtot = -xchg;
		
	      } else if (xchg < -xboxl) {
		fneg=1;
		flag = 1;
		if (-(xchg + xboxl) > xtot)
		  xtot = -(xchg + xboxl);
		
	      }
	    }
	  }
	  if (flag == 1) {
	    
	    /*Put back inside the box, extended out*/
	    traj[0] += 2.0 * xtot;
	    if(fneg>0 && fpos>0){
	      /*For a large complex, test if it could be pushed back out the other side*/
	      cout <<"xtot: "<<xtot<<" but also extends in other direction. negwall "<<negwall<<" poswall: "<<poswall<<endl;
	      recheck=1;
	    }else if(fneg>0){
	      /*Also need to check that update will not push you out the other side*/
	      //xtot is positive. 
	      if(poswall<2.0*xtot){
		recheck=1;
		cout <<" Will push out the other side. negwall "<<negwall<<" poswall: "<<poswall<<" 2*xtot: "<<2.0*xtot<<endl;
	      }
	      
	    }else if(fpos>0){
	      //fpos>0, which means xtot is negative.
	      if(negwall<-2.0*xtot){
		recheck=1;
		cout <<" Will push out the other side. negwall "<<negwall<<" poswall: "<<poswall<<" xtot: "<<xtot<<endl;
	      }
	      
	    }
	    
	    
	  }
	}
	if (flagy >0) {
	  flag = 0;
	  ytot = 0;
	  fneg=0;
	  fpos=0;
	  poswall=plist.yboxl;
	  negwall=plist.yboxl;
	  
	  x0 = ind_com[k].xcom;
	  y0 = ind_com[k].ycom;
	  z0 = ind_com[k].zcom;
	  row[0] = M[3];
	  row[1] = M[4];
	  row[2] = M[5];
	  
	  /*these need to be what current positions
	    due to translation and rotation are*/
	  for (i = 0; i < s1; i++) {
	    mp = ind_com[k].plist[i];
	    
	    /*measure each protein COM to y plane*/
	    dx = bases[mp].xcom - x0;
	    dy = bases[mp].ycom - y0;
	    dz = bases[mp].zcom - z0;
	    /*only need y component*/
	    //rotate(v, M, v2);
	    dyrot = row[0] * dx + row[1] * dy + row[2] * dz;
	    curry = y0 + traj[1] + dyrot;
	    if(yboxl/2.0-curry<poswall)
	      poswall=yboxl/2.0-curry;//shortest distance in y from the +wall
	    if(curry+yboxl/2.0<negwall)
	      negwall=curry+yboxl/2.0;//shortest distance in y from the -wall. >0 inside the box, <0 outside the box.
	    
	    ychg = curry - yboxl / 2.0;
	    if (ychg > 0) {
	      fpos=1;
	      flag = 1;
	      if (-ychg < ytot)
		ytot = -ychg;
	      
	    } else if (ychg < -yboxl) {
	      flag = 1;
	      fneg=1;
	      if (-(ychg + yboxl) > ytot)
		ytot = -(ychg + yboxl);
	    }
	    /*measure each interface to y plane*/
	    for (j = 0; j < bases[mp].ninterface; j++) {
	      dx = bases[mp].x[j] - x0;
	      dy = bases[mp].y[j] - y0;
	      dz = bases[mp].z[j] - z0;
	      /*only need y component*/
	      //rotate(v, M, v2);
	      dyrot = row[0] * dx + row[1] * dy + row[2] * dz;
	      curry = y0 + traj[1] + dyrot;
	      if(yboxl/2.0-curry<poswall)
		poswall=yboxl/2.0-curry;//shortest distance in y from the +wall
	      if(curry+yboxl/2.0<negwall)
		negwall=curry+yboxl/2.0;//shortest distance in y from the -wall. >0 inside the box, <0 outside the box.
	      
	      ychg = curry - yboxl / 2.0;
	      if (ychg > 0) {
		fpos=1;
		flag = 1;
		if (-ychg < ytot)
		  ytot = -ychg;
		
	      } else if (ychg < -yboxl) {
		flag = 1;
		fneg=1;
		if (-(ychg + yboxl) > ytot)
		  ytot = -(ychg + yboxl);
	      }
	    }
	  }
	  if (flag == 1) {
	    
	    /*Put back inside the box*/
	    traj[1] += 2.0 * ytot;
	    if(fneg>0 && fpos>0){
	      /*For a large complex, test if it could be pushed back out the other side*/
	      cout <<"ytot: "<<ytot<<" but also extends in other direction. negwall "<<negwall<<" poswall: "<<poswall<<endl;
	      recheck=1;
	    }else if(fneg>0){
	      /*Also need to check that update will not push you out the other side*/
	      //ytot is positive. 
	      if(poswall<2.0*ytot){
		recheck=1;
		cout <<" Will push out the other side. negwall "<<negwall<<" poswall: "<<poswall<<" ytot: "<<ytot<<endl;
	      }
	      
	    }else if(fpos>0){
	      //fpos>0, which means ytot is negative.
	      if(negwall<-2.0*ytot){
		recheck=1;
		cout <<" Will push out the other side. negwall "<<negwall<<" poswall: "<<poswall<<" ytot: "<<ytot<<endl;
	      }
	      
	    }
		
	  }
	}
	if (flagz >0) {
	  flag = 0;
	  ztot = 0;
	  fneg=0;
	  fpos=0;
	  poswall=plist.yboxl;
	  negwall=plist.yboxl;
	  
	  x0 = ind_com[k].xcom;
	  y0 = ind_com[k].ycom;
	  z0 = ind_com[k].zcom;
	  row[0] = M[6];
	  row[1] = M[7];
	  row[2] = M[8];
	  
	  /*these need to be what current positions
	    due to translation and rotation are*/
	  for (i = 0; i < s1; i++) {
	    mp = ind_com[k].plist[i];
	    
	    /*measure each protein COM to z plane*/
	    dx = bases[mp].xcom - x0;
	    dy = bases[mp].ycom - y0;
	    dz = bases[mp].zcom - z0;
	    
	    dzrot = row[0] * dx + row[1] * dy + row[2] * dz;
	    currz = z0 + traj[2] + dzrot;
	    if(zboxl/2.0-currz<poswall)
	      poswall=zboxl/2.0-currz;//shortest distance in z from the +wall
	    if(currz+zboxl/2.0<negwall)
	      negwall=currz+zboxl/2.0;//shortest distance in z from the -wall. >0 inside the box, <0 outside the box.
	    
	    
	    zchg = currz - zboxl / 2.0;
	    if (zchg > 0) {
	      flag = 1;
	      fpos=1;
	      if (-zchg < ztot)
		ztot = -zchg;
	      
	    } else if (zchg < -zboxl) {
	      flag = 1;
	      fneg=1;
	      if (-(zchg + zboxl) > ztot)
		ztot = -(zchg + zboxl);
	    }
	    /*measure each interface to z plane*/
	    for (j = 0; j < bases[mp].ninterface; j++) {
	      dx = bases[mp].x[j] - x0;
	      dy = bases[mp].y[j] - y0;
	      dz = bases[mp].z[j] - z0;
	      
	      dzrot = row[0] * dx + row[1] * dy + row[2] * dz;
	      currz = z0 + traj[2] + dzrot;
	      if(zboxl/2.0-currz<poswall)
		poswall=zboxl/2.0-currz;//shortest distance in z from the +wall
	      if(currz+zboxl/2.0<negwall)
		negwall=currz+zboxl/2.0;//shortest distance in z from the -wall. >0 inside the box, <0 outside the box.
	      
	      
	      zchg = currz - zboxl / 2.0;
	      if (zchg > 0) {
		flag = 1;
		fpos=1;
		if (-zchg < ztot)
		  ztot = -zchg;
		
	      } else if (zchg < -zboxl) {
		flag = 1;
		fneg=1;
		if (-(zchg + zboxl) > ztot)
		  ztot = -(zchg + zboxl);
	      }
	    }
	  }
	  if (flag == 1) {
	    
	    /*Put back inside the box*/
	    traj[2] += 2.0 * ztot;
	    if(fneg>0 && fpos>0){
	      /*For a large complex, test if it could be pushed back out the other side*/
	      cout <<"ztot: "<<ztot<<" but also extends in other direction. negwall "<<negwall<<" poswall: "<<poswall<<endl;
	      recheck=1;
	    }else if(fneg>0){
	      /*Also need to check that update will not push zou out the other side*/
	      //ztot is positive. 
	      if(poswall<2.0*ztot){
		recheck=1;
		cout <<" Will push out the other side. negwall "<<negwall<<" poswall: "<<poswall<<" ztot: "<<ztot<<endl;
	      }
	    }else if(fpos>0){
	      //fpos>0, which means ztot is negative.
	      if(negwall<-2.0*ztot){
		recheck=1;
		cout <<" Will push out the other side. negwall "<<negwall<<" poswall: "<<poswall<<" ztot: "<<ztot<<endl;
	      }
	      
	    }
	

	  }//updating traj to reflect
	}
	if(recheck>0){
	  /*Test that new coordinates have not pushed you out of the box for a very large complex, if so, resample rotation matrix.*/
	  cout <<"RECHECK THAT COMPLEX DOES NOT SPAN BOX IN SUBROUTINE REFLECT TRAJ COMPLEX RAD ROT. Complex: "<<k<<" size:" <<s1<<endl;
	  reflect_traj_check_span( xtot,  ytot, ztot, p1,  bases, ind_com, plist, traj, trajR, M) ;
	  
	}
}
