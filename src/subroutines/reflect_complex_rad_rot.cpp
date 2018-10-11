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

void reflect_complex_rad_rot(int p1, Fullmol *bases, Complex *ind_com, Parms plist) {
	int i, j;
	int k = bases[p1].mycomplex;
	int s1 = ind_com[k].mysize;
	double xchg;
	double ychg;
	double zchg;
	double xtot = 0.0;
	double ytot = 0.0;
	double ztot = 0.0;
	int flag = 0;
	double rad = ind_com[k].radR;
	int flagz = 0;
	int flagy = 0;
	int flagx = 0;
	double xtot0, ytot0, ztot0;
	double tol = 1E-11;
	double currx, curry, currz;
	double poswall, negwall;
	double xboxl=plist.xboxl;
	double yboxl=plist.yboxl;
	double zboxl=plist.zboxl;
	xchg = ind_com[k].xcom - xboxl / 2.0;
	if ((xchg + rad) > 0) 
		flagx = 1;
	if ((xchg - rad) < -xboxl) 
		flagx += 1;
	ychg = ind_com[k].ycom - yboxl / 2.0;
	if ((ychg + rad) > 0) 
		flagy = 1;
	 if ((ychg - rad) < -yboxl) 
		flagy += 1;
	 zchg = ind_com[k].zcom - zboxl / 2.0;
	if ((zchg + rad) > 0) 
		flagz = 1;
	 if ((zchg - rad) < -zboxl) 
		flagz += 1;
	
	xtot=0;
	ytot=0;
	ztot=0;
	int mp;
	int fpos, fneg;
	if (flagz >0) {
	  poswall=zboxl;
	  negwall=zboxl;
	  fneg=0;
	  fpos=0;
	  flag = 0;
	  ztot = 0;
	  for (i = 0; i < s1; i++) {
	    mp = ind_com[k].plist[i];
	    
	    /*also measure each protein COM to z-plane*/
	    zchg = bases[mp].zcom - zboxl / 2.0;
	    currz=bases[mp].zcom;
	    if(zboxl/2.0-currz<poswall)
	      poswall=zboxl/2.0-currz;//shortest distance in z from the +wall
	    if(currz+zboxl/2.0<negwall)
	      negwall=currz+zboxl/2.0;//shortest distance in z from the -wall. >0 inside the box, <0 outside the box.
	    
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
	      zchg = bases[mp].z[j] - zboxl / 2.0;
	      currz=bases[mp].z[j];
	      if(zboxl/2.0-currz<poswall)
		poswall=zboxl/2.0-currz;//shortest distance in z from the +wall
	      if(currz+zboxl/2.0<negwall)
		negwall=currz+zboxl/2.0;//shortest distance in z from the -wall. >0 inside the box, <0 outside the box.
	      
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
	    
	    if(fneg>0 && fpos>0){
	      //extends out both the front and back.
	      cout <<"IN REFLECT COMPLEX RAD ROT, EXTENDS . ALREADY UPDATED POSITIONS. EXITING..."<<endl;
	      cout <<"ztot: "<<ztot<<" but also extends in other direction. negwall "<<negwall<<" poswall: "<<poswall<<endl;
	      exit(1);
	    }else if(fneg>0){
	      /*Also need to check that update will not push zou out the other side*/
	      //ztot is positive. 
	      if(poswall<2.0*ztot){
		cout <<" IN REFLECT COMPLEX RAD ROT: Will push out the other side. negwall "<<negwall<<" poswall: "<<poswall<<" ztot: "<<ztot<<endl;
		if(poswall<ztot){
		  cout <<" PROBLEM...IN REFLECT COMPLEX RAD ROT, EXTENDS BOTH OUT FRONT AND BACK.."<<endl;
		  exit(1);
		}else{
		  /*Fit back inside boz, without bouncing off (just at boz edge)*/
		  ind_com[k].zcom += ztot;
		  for (i = 0; i < s1; i++) {
		    mp = ind_com[k].plist[i];
		    bases[mp].zcom +=  ztot;
		    //update interface coords
		    for (j = 0; j < bases[mp].ninterface; j++) {
		      bases[mp].z[j] +=  ztot;
		    }
		  }
		  
		}

	      }else{
		cout <<"IN REFLECT COMPLEX RAD ROT: put back in the box, z: "<<ztot<<endl;
		//just update positions.
		/*Put back inside the box*/
		ind_com[k].zcom += 2.0 * ztot;
		
		//update protein COM
		
		for (i = 0; i < s1; i++) {
		  mp = ind_com[k].plist[i];
		  bases[mp].zcom += 2.0 * ztot;
		  //update interface coords
		  for (j = 0; j < bases[mp].ninterface; j++) {
		    bases[mp].z[j] += 2.0 * ztot;
		  }
		}
		
	      }
	      
	    }else{
	      //fpos>0, which means ztot is negative.
	      if(negwall<-2.0*ztot){
		
		cout <<" IN REFLECT COMPLEX RAD ROT: Will push out the other side. negwall "<<negwall<<" poswall: "<<poswall<<" ztot: "<<ztot<<endl;
		if(negwall<-ztot){
		  cout <<" PROBLEM...IN REFLECT COMPLEX RAD ROT, EXTENDS BOTH OUT FRONT AND BACK.."<<endl;
		  exit(1);
		}else{
		  /*Fit back inside boz, without bouncing off (just at boz edge)*/
		  ind_com[k].zcom += ztot;
		  for (i = 0; i < s1; i++) {
		    mp = ind_com[k].plist[i];
		    bases[mp].zcom +=  ztot;
		    //update interface coords
		    for (j = 0; j < bases[mp].ninterface; j++) {
		      bases[mp].z[j] +=  ztot;
		    }
		  }
		  
		}
	      }else{
		cout <<"IN REFLECT COMPLEX RAD ROT: put back in the box, z: "<<ztot<<endl;
		//just update positions.
		/*Put back inside the box*/
		ind_com[k].zcom += 2.0 * ztot;
		//update protein COM
		for (i = 0; i < s1; i++) {
		  mp = ind_com[k].plist[i];
		  bases[mp].zcom += 2.0 * ztot;
		  //update interface coords
		  for (j = 0; j < bases[mp].ninterface; j++) {
		    bases[mp].z[j] += 2.0 * ztot;
		  }
		}
		
	      }
	      
	    }
	    	    
	  }//update traj for reflection
	  
	}
	if (flagy >0) {
	  
	  flag = 0;
	  ytot = 0;
	  fneg=0;
	  fpos=0;
	  poswall=yboxl;
	  negwall=yboxl;
	  for (i = 0; i < s1; i++) {
	    mp = ind_com[k].plist[i];
	    /*measure each protein COM to y plane*/
	    curry=bases[mp].ycom;
	    if(yboxl/2.0-curry<poswall)
	      poswall=yboxl/2.0-curry;//shortest distance in y from the +wall
	    if(curry+yboxl/2.0<negwall)
	      negwall=curry+yboxl/2.0;//shortest distance in y from the -wall. >0 inside the box, <0 outside the box.
	    
	    ychg = bases[mp].ycom - yboxl / 2.0;
	    
	    if (ychg > 0) {
	      flag = 1;
	      fpos=1;
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
	      curry=bases[mp].y[j];
	      if(yboxl/2.0-curry<poswall)
		poswall=yboxl/2.0-curry;//shortest distance in y from the +wall
	      if(curry+yboxl/2.0<negwall)
		negwall=curry+yboxl/2.0;//shortest distance in y from the -wall. >0 inside the box, <0 outside the box.
	      
	      ychg = bases[mp].y[j] - yboxl / 2.0;
	      
	      if (ychg > 0) {
		flag = 1;
		fpos=1;
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
	    
	    if(fneg>0 && fpos>0){
	      //extends out both the front and back.
	      cout <<"IN REFLECT COMPLEX RAD ROT, EXTENDS . ALREADY UPDATED POSITIONS, EXITING... "<<endl;
	      cout <<"ytot: "<<ytot<<" but also extends in other direction. negwall "<<negwall<<" poswall: "<<poswall<<endl;
	      exit(1);
	    }else if(fneg>0){
	      /*Also need to check that update will not push you out the other side*/
	      //ytot is positive. 
	      if(poswall<2.0*ytot){
		cout <<" IN REFLECT COMPLEX RAD ROT: Will push out the other siye. negwall "<<negwall<<" poswall: "<<poswall<<" ytot: "<<ytot<<endl;
		if(poswall<ytot){
		  cout <<" PROBLEM...IN REFLECT COMPLEX RAD ROT, EXTENDS BOTH OUT FRONT AND BACK.."<<endl;
		  exit(1);
		}else{
		  /*Fit back inside boy, without bouncing off (just at boy edge)*/
		  ind_com[k].ycom += ytot;
		  for (i = 0; i < s1; i++) {
		    mp = ind_com[k].plist[i];
		    bases[mp].ycom +=  ytot;
		    //update interface coords
		    for (j = 0; j < bases[mp].ninterface; j++) {
		      bases[mp].y[j] +=  ytot;
		    }
		  }
		  
		}

	      }else{
		cout <<"IN REFLECT COMPLEX RAD ROT: put back in the box, y: "<<ytot<<endl;
		//just update positions.
		/*Put back inside the box*/
		ind_com[k].ycom += 2.0 * ytot;
		
		//update protein COM
		
		for (i = 0; i < s1; i++) {
		  mp = ind_com[k].plist[i];
		  bases[mp].ycom += 2.0 * ytot;
		  //update interface coords
		  for (j = 0; j < bases[mp].ninterface; j++) {
		    bases[mp].y[j] += 2.0 * ytot;
		  }
		}
		
	      }
	      
	    }else{
	      //fpos>0, which means ytot is negative.
	      if(negwall<-2.0*ytot){
		
		cout <<" IN REFLECT COMPLEX RAD ROT: Will push out the other siye. negwall "<<negwall<<" poswall: "<<poswall<<" ytot: "<<ytot<<endl;
		if(negwall<-ytot){
		  cout <<" PROBLEM...IN REFLECT COMPLEX RAD ROT, EXTENDS BOTH OUT FRONT AND BACK.."<<endl;
		  exit(1);
		}else{
		  /*Fit back inside boy, without bouncing off (just at boy edge)*/
		  ind_com[k].ycom += ytot;
		  for (i = 0; i < s1; i++) {
		    mp = ind_com[k].plist[i];
		    bases[mp].ycom +=  ytot;
		    //update interface coords
		    for (j = 0; j < bases[mp].ninterface; j++) {
		      bases[mp].y[j] +=  ytot;
		    }
		  }
		  
		}

	      }else{
		cout <<"IN REFLECT COMPLEX RAD ROT: put back in the box, y: "<<ytot<<endl;
		//just update positions.
		/*Put back inside the box*/
		ind_com[k].ycom += 2.0 * ytot;
		//update protein COM
		for (i = 0; i < s1; i++) {
		  mp = ind_com[k].plist[i];
		  bases[mp].ycom += 2.0 * ytot;
		  //update interface coords
		  for (j = 0; j < bases[mp].ninterface; j++) {
		    bases[mp].y[j] += 2.0 * ytot;
		  }
		}
		
	      }
	      
	    }
	    	    
	  }//update traj for reflection
	}
	if (flagx >0) {

	  fneg=0;
	  fpos=0;
	  negwall=xboxl;
	  poswall=xboxl;
	  flag = 0;
	  xtot = 0;
	  for (i = 0; i < s1; i++) {
	    mp = ind_com[k].plist[i];
	    
	    /*measure each protein COM to x-plane*/
	    currx=bases[mp].xcom;
	    if(xboxl/2.0-currx<poswall)
	      poswall=xboxl/2.0-currx;//shortest distance in x from the +wall
	    if(currx+xboxl/2.0<negwall)
	      negwall=currx+xboxl/2.0;//shortest distance in x from the -wall. >0 inside the box, <0 outside the box.
	    
	    xchg = bases[mp].xcom - xboxl / 2.0;
	    
	    if (xchg > 0) {
	      flag = 1;
	      fpos=1;
	      if (-xchg < xtot)
		xtot = -xchg;
	      
	    } else if (xchg < -xboxl) {
	      fneg=1;
	      flag = 1;
	      if (-(xchg + xboxl) > xtot)
		xtot = -(xchg + xboxl);
	      
	    }
	    /*measure each interface to x plane*/
	    for (j = 0; j < bases[mp].ninterface; j++) {
	      currx=bases[mp].x[j];
	      if(xboxl/2.0-currx<poswall)
		poswall=xboxl/2.0-currx;//shortest distance in x from the +wall
	      if(currx+xboxl/2.0<negwall)
		negwall=currx+xboxl/2.0;//shortest distance in x from the -wall. >0 inside the box, <0 outside the box.
	      
	      xchg = bases[mp].x[j] - xboxl / 2.0;
	      
	      if (xchg > 0) {
		flag = 1;
		fpos=1;
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
	    /*Move inside box, but test for box span.*/
	    if(fneg>0 && fpos>0){
	      //extends out both the front and back.
	      cout <<"IN REFLECT COMPLEX RAD ROT, EXTENDS . ALREADY MOVED WITHOUT TRAJ, EXITING.... "<<endl;
	      cout <<"xtot: "<<xtot<<" but also extends in other direction. negwall "<<negwall<<" poswall: "<<poswall<<endl;
	      exit(1);
	    }else if(fneg>0){
	      /*Also need to check that update will not push xou out the other side*/
	      //xtot is positive. 
	      if(poswall<2.0*xtot){
		cout <<" IN REFLECT COMPLEX RAD ROT: Will push out the other sixe. negwall "<<negwall<<" poswall: "<<poswall<<" xtot: "<<xtot<<endl;
		if(poswall<xtot){
		  cout <<" PROBLEM...IN REFLECT COMPLEX RAD ROT, EXTENDS BOTH OUT FRONT AND BACK.."<<endl;
		  exit(1);
		}else{
		  /*Fit back inside box, without bouncing off (just at box edge)*/
		  ind_com[k].xcom += xtot;
		  for (i = 0; i < s1; i++) {
		    mp = ind_com[k].plist[i];
		    bases[mp].xcom +=  xtot;
		    //update interface coords
		    for (j = 0; j < bases[mp].ninterface; j++) {
		      bases[mp].x[j] +=  xtot;
		    }
		  }
		  
		}
		

	      }else{
		cout <<"IN REFLECT COMPLEX RAD ROT: put back in the box, x: "<<xtot<<endl;
		//just update positions.
		/*Put back inside the box*/
		ind_com[k].xcom += 2.0 * xtot;
		
		//update protein COM
		
		for (i = 0; i < s1; i++) {
		  mp = ind_com[k].plist[i];
		  bases[mp].xcom += 2.0 * xtot;
		  //update interface coords
		  for (j = 0; j < bases[mp].ninterface; j++) {
		    bases[mp].x[j] += 2.0 * xtot;
		  }
		}
		
	      }
	      
	    }else{
	      //fpos>0, which means xtot is negative.
	      if(negwall<-2.0*xtot){
		
		cout <<" IN REFLECT COMPLEX RAD ROT: Will push out the other sixe. negwall "<<negwall<<" poswall: "<<poswall<<" xtot: "<<xtot<<endl;
		if(negwall<-xtot){
		  cout <<" PROBLEM...IN REFLECT COMPLEX RAD ROT, EXTENDS BOTH OUT FRONT AND BACK.."<<endl;
		  exit(1);
		}else{
		  /*Fit back inside box, without bouncing off (just at box edge)*/
		  ind_com[k].xcom += xtot;
		  for (i = 0; i < s1; i++) {
		    mp = ind_com[k].plist[i];
		    bases[mp].xcom +=  xtot;
		    //update interface coords
		    for (j = 0; j < bases[mp].ninterface; j++) {
		      bases[mp].x[j] +=  xtot;
		    }
		  }
		  
		}

	      }else{
		cout <<"IN REFLECT COMPLEX RAD ROT: put back in the box, x: "<<xtot<<endl;
		//just update positions.
		/*Put back inside the box*/
		ind_com[k].xcom += 2.0 * xtot;
		//update protein COM
		for (i = 0; i < s1; i++) {
		  mp = ind_com[k].plist[i];
		  bases[mp].xcom += 2.0 * xtot;
		  //update interface coords
		  for (j = 0; j < bases[mp].ninterface; j++) {
		    bases[mp].x[j] += 2.0 * xtot;
		  }
		}
		
	      }
	      
	    }
	    	    
	  }//update traj for reflection

		
	
		
	}
	

}
