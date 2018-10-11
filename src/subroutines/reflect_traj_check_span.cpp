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

void reflect_traj_check_span(double xtot, double ytot, double ztot, int p1, Fullmol *bases, Complex *ind_com, Parms plist, double *traj, double *trajR, double *M) {
  int i, j;
  int k = bases[p1].mycomplex;
  int s1 = ind_com[k].mysize;
  double xboxl=plist.xboxl;
  double yboxl=plist.yboxl;
  double zboxl=plist.zboxl;
  double deltat=plist.deltat;
  double xchg;
  double ychg;
  double zchg;
  
  double rad = ind_com[k].radR;

  double currx;
  double curry;
  double currz;

  
  double tx, ty, tz;
  int mp;
  int recheck=0;
  double dx, dy, dz;
  double x0, y0, z0;
  double dzrot, dxrot, dyrot;
  double *v=new double[3];
  double *vr=new double[3];
  int flag=1;
  int flagdim=0;
  int maxit=50;
  int it=0;
  int fposx, fposy, fposz;
  int fnegx, fnegy, fnegz;
  double poswallx, poswally, poswallz;
  double negwallx, negwally, negwallz;
  
  double tol=1E-14;
  int fail=0;
  while(it<maxit && flag==1){
    flag=0;//without double span, this will stay 0
    flagdim=0;
    fail=0;
    x0 = ind_com[k].xcom;
    y0 = ind_com[k].ycom;
    z0 = ind_com[k].zcom;
    poswallx=plist.xboxl;
    negwallx=plist.xboxl;
    fposx=0;
    fnegx=0;
    poswally=plist.yboxl;
    negwally=plist.yboxl;
    fposy=0;
    fnegy=0;
    poswallz=plist.zboxl;
    negwallz=plist.zboxl;
    fposz=0;
    fnegz=0;

    xtot=0;
    ytot=0;
    ztot=0;
    
    /*these need to be what current positions
      due to translation and rotation are*/
    for (i = 0; i < s1; i++) {
      mp = ind_com[k].plist[i];
      /*measure protein COM to plane*/
      v[0] = bases[mp].xcom - x0;
      v[1] = bases[mp].ycom - y0;
      v[2] = bases[mp].zcom - z0;
      
      rotate(v, M, vr);
      currx = x0 + traj[0] + vr[0];
      curry=y0+traj[1]+vr[1];
      currz=z0+traj[2]+vr[2];
		
      if(xboxl/2.0-currx<poswallx)
	poswallx=xboxl/2.0-currx;//shortest distance in x from the +wall
      if(currx+xboxl/2.0<negwallx)
	negwallx=currx+xboxl/2.0;//shortest distance in x from the -wall. >0 inside the box, <0 outside the box.
      
      xchg = currx - xboxl / 2.0;
      if (xchg > 0) {
	fposx=1;
	flagdim = 1;
	if (-xchg < xtot)
	  xtot = -xchg;  
      } else if (xchg < -xboxl) {
	fnegx=1;
	flagdim = 1;
	if (-(xchg + xboxl) > xtot)
	  xtot = -(xchg + xboxl);
      }
      
      if(yboxl/2.0-curry<poswally)
	poswally=yboxl/2.0-curry;//shortest distance in y from the +wall
      if(curry+yboxl/2.0<negwally)
	negwally=curry+yboxl/2.0;//shortest distance in y from the -wall. >0 inside the boy, <0 outside the boy.
      
      ychg = curry - yboxl / 2.0;
      if (ychg > 0) {
	fposy=1;
	flagdim = 1;
	if (-ychg < ytot)
	  ytot = -ychg;
	
      } else if (ychg < -yboxl) {
	flagdim = 1;
	fnegy=1;
	if (-(ychg + yboxl) > ytot)
	  ytot = -(ychg + yboxl);
      }
      
      if(zboxl/2.0-currz<poswallz)
	poswallz=zboxl/2.0-currz;//shortest distance in z from the +wall
      if(currz+zboxl/2.0<negwallz)
	negwallz=currz+zboxl/2.0;//shortest distance in z from the -wall. >0 inside the boz, <0 outside the boz.
      
      zchg = currz - zboxl / 2.0;
      if (zchg > 0) {
	fposz=1;
	flagdim = 1;
	if (-zchg < ztot)
	  ztot = -zchg;
	
      } else if (zchg < -zboxl) {
	flagdim = 1;
	fnegz=1;
	if (-(zchg + zboxl) > ztot)
	  ztot = -(zchg + zboxl);
      }
      
      /*measure each interface to z plane*/
      for (j = 0; j < bases[mp].ninterface; j++) {
	v[0] = bases[mp].x[j] - x0;
	v[1] = bases[mp].y[j] - y0;
	v[2] = bases[mp].z[j] - z0;

	rotate(v, M, vr);
	currx = x0 + traj[0] + vr[0];
	curry=y0+traj[1]+vr[1];
	currz=z0+traj[2]+vr[2];
		
	if(xboxl/2.0-currx<poswallx)
	  poswallx=xboxl/2.0-currx;//shortest distance in x from the +wall
	if(currx+xboxl/2.0<negwallx)
	  negwallx=currx+xboxl/2.0;//shortest distance in x from the -wall. >0 inside the box, <0 outside the box.
	      
	xchg = currx - xboxl / 2.0;
	if (xchg > 0) {
	  fposx=1;
	  flagdim = 1;
	  if (-xchg < xtot)
	    xtot = -xchg;  
	} else if (xchg < -xboxl) {
	  fnegx=1;
	  flagdim = 1;
	  if (-(xchg + xboxl) > xtot)
	    xtot = -(xchg + xboxl);
	}

	if(yboxl/2.0-curry<poswally)
	  poswally=yboxl/2.0-curry;//shortest distance in y from the +wall
	if(curry+yboxl/2.0<negwally)
	  negwally=curry+yboxl/2.0;//shortest distance in y from the -wall. >0 inside the boy, <0 outside the boy.
	
	ychg = curry - yboxl / 2.0;
	if (ychg > 0) {
	  fposy=1;
	  flagdim = 1;
	  if (-ychg < ytot)
	    ytot = -ychg;
	  
	} else if (ychg < -yboxl) {
	  flagdim = 1;
	  fnegy=1;
	  if (-(ychg + yboxl) > ytot)
	    ytot = -(ychg + yboxl);
	}

	if(zboxl/2.0-currz<poswallz)
	  poswallz=zboxl/2.0-currz;//shortest distance in z from the +wall
	if(currz+zboxl/2.0<negwallz)
	  negwallz=currz+zboxl/2.0;//shortest distance in z from the -wall. >0 inside the boz, <0 outside the boz.
	
	zchg = currz - zboxl / 2.0;
	if (zchg > 0) {
	  fposz=1;
	  flagdim = 1;
	  if (-zchg < ztot)
	    ztot = -zchg;
	  
	} else if (zchg < -zboxl) {
	  flagdim = 1;
	  fnegz=1;
	  if (-(zchg + zboxl) > ztot)
	    ztot = -(zchg + zboxl);
	}

      }//loop over interfaces
    }//loop over proteins in complex.
    
    if (flagdim==1) {
      /*Put back inside the box, extended out*/
      traj[0] += 2.0 * xtot;
      traj[1] += 2.0 * ytot;
      traj[2] += 2.0 * ztot;
      if(fnegx>0 && fposx>0){
	/*For a large complex, test if it could be pushed back out the other side*/
	cout <<"xtot: "<<xtot<<" but also extends in other direction. negwall "<<negwallx<<" poswall: "<<poswallx<<endl;
	cout <<"Pushed back out of box X, resample M and traj "<<endl;
	fail=1;
	
      }else if(fnegx>0){
	/*Also need to check that update will not push you out the other side*/
	//xtot is positive. 
	if(poswallx<2.0*xtot){
	  cout <<"Pushed back out of box X, resample M and traj "<<endl;
	  fail=1;
	  cout <<" Will push out the other side. negwall "<<negwallx<<" poswall: "<<poswallx<<" 2*xtot: "<<2.0*xtot<<endl;
	}
	
      }else if(fposx>0){
	//fpos>0, which means xtot is negative.
	if(negwallx<-2.0*xtot){
	  cout <<"Pushed back out of box X, resample M and traj "<<endl;
	  fail=1;
	  cout <<" Will push out the other side. negwall "<<negwallx<<" poswall: "<<poswallx<<" xtot: "<<xtot<<endl;
	}
	
      }
      if(fnegy>0 && fposy>0){
	/*For a large complex, test if it could be pushed back out the other side*/
	cout <<"ytot: "<<ytot<<" but also extends in other direction. negwall "<<negwally<<" poswall: "<<poswally<<endl;
	cout <<"Pushed back out of box Y, resample M and traj "<<endl;
	fail=1;
      }else if(fnegy>0){
	/*Also need to check that update will not push you out the other side*/
	//xtot is positive. 
	if(poswally<2.0*ytot){
	  cout <<"Pushed back out of box Y, resample M and traj "<<endl;
	  fail=1;
	  cout <<" Will push out the other side. negwall "<<negwally<<" poswall: "<<poswally<<" 2*ytot: "<<2.0*ytot<<endl;
	}
	
      }else if(fposy>0){
	//fpos>0, which means xtot is negative.
	if(negwally<-2.0*ytot){
	  cout <<"Pushed back out of box, Y resample M and traj "<<endl;
	  fail=1;
	  cout <<" Will push out the other side. negwall "<<negwally<<" poswall: "<<poswally<<" ytot: "<<ytot<<endl;
	}
	
      }
      if(fnegz>0 && fposz>0){
	/*For a large complex, test if it could be pushed back out the other side*/
	cout <<"ztot: "<<ztot<<" but also extends in other direction. negwall "<<negwallz<<" poswall: "<<poswallz<<endl;
	cout <<"Pushed back out of box Z, resample M and traj "<<endl;
	fail=1;
      }else if(fnegz>0){
	/*Also need to check that update will not push zou out the other side*/
	//xtot is positive. 
	if(poswallz<2.0*ztot){
	  cout <<"Pushed back out of box Z, resample M and traj "<<endl;
	  fail=1;
	  cout <<" Will push out the other side. negwall "<<negwallz<<" poswall: "<<poswallz<<" 2*ztot: "<<2.0*ztot<<endl;
	}
	
      }else if(fposz>0){
	//fpos>0, which means xtot is negative.
	if(negwallz<-2.0*ztot){
	  cout <<"Pushed back out of box, Z resample M and traj "<<endl;
	  fail=1;
	}
      }
    }//recheck span
    if(fail==1){
      /*Resample, extends in x, y, and/or z*/
      traj[0] = sqrt(2.0 * deltat * ind_com[k].Dx) * GaussV();
      traj[1] = sqrt(2.0 * deltat * ind_com[k].Dy) * GaussV();
      traj[2] = sqrt(2.0 * deltat * ind_com[k].Dz) * GaussV();
      trajR[0] = sqrt(2.0 * deltat * ind_com[k].Drx) * GaussV();
      trajR[1] = sqrt(2.0 * deltat * ind_com[k].Dry) * GaussV();
      trajR[2] = sqrt(2.0 * deltat * ind_com[k].Drz) * GaussV();
      
      rotationEuler(trajR[0], trajR[1], trajR[2], M);
      reflect_traj_complex_rad_rot_nocheck(p1, bases, ind_com, plist.xboxl, plist.yboxl, plist.zboxl, traj, M);
      it++;
      flag=1;//will need to recheck after resampling traj and trajR
    }
  }//loop over iterations and flag condition
  
  cout <<"ITERATIONS TO CONVERGE POSITION WITHIN BOX: "<<it<<" flag at end: 0=success: "<<flag<<endl;
  if(flag==1){
    cout <<"WARNING: DID NOT CONVERGE POSITION, NEW POS: "<<endl;
    
    for (i = 0; i < s1; i++) {
      mp = ind_com[k].plist[i];
      v[0] = bases[mp].xcom - x0;
      v[1] = bases[mp].ycom - y0;
      v[2] = bases[mp].zcom - z0;
      double dx=traj[0];
      double dy=traj[1];
      double dz=traj[2];
      rotate(v, M, vr);
      /*first would make xcom=x0+vr, then would also add dx */
      cout <<"i: "<<i<<" P: "<<mp<<" com:" <<x0+dx+vr[0]<<' '<<y0+dy+vr[1]<<' '<<z0+dz+vr[2]<<endl;
      
      //update interface coords
      for (j = 0; j < bases[mp].ninterface; j++) {
	v[0] = bases[mp].x[j] - x0;
	v[1] = bases[mp].y[j] - y0;
	v[2] = bases[mp].z[j] - z0;
	rotate(v, M, vr);
	/*first would make xcom=x0+vr, then would also add dx */
	cout <<x0+dx+vr[0]<<' '<<y0+dy+vr[1]<<' '<<z0+dz+vr[2]<<endl;
      }
    }
  }//DID NOT CONVERGE
  delete[] v;
  delete[] vr;
  
}
