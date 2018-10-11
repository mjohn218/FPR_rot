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

int check_box_span(int p1, int p2, Fullmol *bases, Complex *ind_com, Parms &plist){
  /*Associating proteins have been moved to contact. Before assigning them to the same complex,
    test to see if the complex is too big to fit in the box.
   */
  int i, j;
  double xboxl=plist.xboxl;
  double yboxl=plist.yboxl;
  double zboxl=plist.zboxl;
  int k = bases[p1].mycomplex;
  int s1 = ind_com[k].mysize;
  int k2=bases[p2].mycomplex;
  int s2=ind_com[k2].mysize;
  double xchg;
  double ychg;
  double zchg;
  double currx, curry, currz;
  double xtot = 0.0;
  double ytot = 0.0;
  double ztot = 0.0;
  int flag = 0;
  double rad = ind_com[k].radR;
  double rad2 = ind_com[k2].radR;
  int flagz = 0;
  int flagy = 0;
  int flagx = 0;
  double xtot2, ytot2, ztot2;
  double tol = 1E-11;
  int cancel=0;
  //Complex 1, x

  int fpos=0;
  int fneg=0;
  double negwall, poswall;
  if(rad+rad2>xboxl/2.0)flagx=2;
  if(rad+rad2>yboxl/2.0)flagy=2;
  if(rad+rad2>zboxl/2.0)flagz=2;

  int mp;
  int flagpos, flagneg;
  if (flagz >1) {
    cout <<"CHECK BOX SPAN in Z. Complex 1 Radius, Complex 2 radius "<<rad<<' '<<rad2<<endl;
    //The approximate size of the complex in z (max size) puts it as outside, now test interface positions.
    flagpos = 0;
    flagneg=0;
    ztot = 0;
    ztot2=0;
    poswall=plist.zboxl;
    negwall=plist.zboxl;
    
    for (i = 0; i < s1; i++) {
	    mp = ind_com[k].plist[i];
	    /*measure each interface to z plane*/
	    for (j = 0; j < bases[mp].ninterface; j++) {
	      currz = bases[mp].z[j];
	      if(zboxl/2.0-currz<poswall)
		poswall=zboxl/2.0-currz;//shortest distance in z from the +wall
	      if(currz+zboxl/2.0<negwall)
		negwall=currz+zboxl/2.0;//shortest distance in z from the -wall. >0 inside the box, <0 outside the box.
	     
	      zchg = bases[mp].z[j] - zboxl / 2.0;
	      if (zchg > 0) {
		flagpos = 1;
		if (-zchg < ztot)
		  ztot = -zchg;
	      }
	      if (zchg < -zboxl) {
		flagneg = 1;
		if (-(zchg + zboxl) > ztot)
		  ztot = -(zchg + zboxl);
	      }
	    }
	  }//finished collecting flagpos and flagneg for complex 1
	  for (i = 0; i < s2; i++) {
	    mp = ind_com[k2].plist[i];
	    /*measure each interface to z plane*/
	    for (j = 0; j < bases[mp].ninterface; j++) {
	      currz = bases[mp].z[j];
	      if(zboxl/2.0-currz<poswall)
		poswall=zboxl/2.0-currz;//shortest distance in z from the +wall
	      if(currz+zboxl/2.0<negwall)
		negwall=currz+zboxl/2.0;//shortest distance in z from the -wall. >0 inside the box, <0 outside the box.
	     
	      zchg = bases[mp].z[j] - zboxl / 2.0;
	      if (zchg > 0) {
		flagpos = 1;
		if (-zchg < ztot2)
		  ztot2 = -zchg;
	      }
	      if (zchg < -zboxl) {
		flagneg = 1;
		if (-(zchg + zboxl) > ztot2)
		  ztot2 = -(zchg + zboxl);
	      }
	    }
	  }//finished collecting flagpos and flagneg, for complex 2.
	  //cout <<"flagneg and flagpos:" <<flagneg<<' '<<flagpos<<endl;
	  //cout <<"poswal and negwall: "<<poswall<<' '<<negwall<<endl;
	  
	  if (flagneg >0 && flagpos>0) {
	    cancel=-1;//current pos spans both
	    cout <<"STICKS OUT BOTH SIZES in Z, CANCEL ASSOCIATION "<<endl;
	    //write_crds_complex(k, ind_com, bases);
	    //write_crds_complex(k2, ind_com, bases);
	  }else if(flagneg>0){
	    //ztot is positive.
	    if(poswall<ztot || poswall<ztot2){
	      cancel=-1;//current pos spans both
	      cout <<"WILL STICK OUT in Z, CANCEL ASSOCIATION "<<endl;
	      //write_crds_complex(k, ind_com, bases);
	      //write_crds_complex(k2, ind_com, bases);
	    }else{
	      //put back in the box. put at edge, rather than bouncing off.
	      if(ztot2>ztot)ztot=ztot2;
	      ind_com[k].zcom +=  ztot;
	      for (i = 0; i < s1; i++) {
		mp = ind_com[k].plist[i];
		bases[mp].zcom +=  ztot;
		  //update interface coords
		  for (j = 0; j < bases[mp].ninterface; j++) {
		    bases[mp].z[j] += ztot;
		  }
	      }
	      ind_com[k2].zcom += ztot;
	      for (i = 0; i < s2; i++) {
		  mp = ind_com[k2].plist[i];
		  bases[mp].zcom +=  ztot;
		  //update interface coords
		  for (j = 0; j < bases[mp].ninterface; j++) {
		    bases[mp].z[j] +=  ztot;
		  }
		}
		
	    }
	  }else if(flagpos>0){
	    //fpos>0, which means ztot is negative.
	    if(negwall<-ztot || negwall<-ztot2){
	      cancel=-1;//current pos spans both
	      cout <<"WILL STICK OUT in Z, CANCEL ASSOCIATION "<<endl;
	      //write_crds_complex(k, ind_com, bases);
	      //write_crds_complex(k2, ind_com, bases);
	    }else{
	      //put back in the box. put at edge, rather than bouncing off.
	      if(ztot2<ztot)ztot=ztot2;
	      ind_com[k].zcom +=  ztot;
	      for (i = 0; i < s1; i++) {
		mp = ind_com[k].plist[i];
		bases[mp].zcom +=  ztot;
		  //update interface coords
		  for (j = 0; j < bases[mp].ninterface; j++) {
		    bases[mp].z[j] += ztot;
		  }
	      }
	      ind_com[k2].zcom += ztot;
	      for (i = 0; i < s2; i++) {
		  mp = ind_com[k2].plist[i];
		  bases[mp].zcom +=  ztot;
		  //update interface coords
		  for (j = 0; j < bases[mp].ninterface; j++) {
		    bases[mp].z[j] +=  ztot;
		  }
		}

	    }
	  }
	}//flagz>1
	  
	if (flagy >1) {
	  cout <<"CHECK BOX SPAN in Y. Complex 1 Radius, Complex 2 radius "<<rad<<' '<<rad2<<endl;
	  poswall=plist.yboxl;
	  negwall=plist.yboxl;
	  
	  flagpos = 0;
	  flagneg=0;
	  ytot = 0;
	  ytot2=0;
	  for (i = 0; i < s1; i++) {
	    mp = ind_com[k].plist[i];
	    /*measure each interface to y plane*/
	    for (j = 0; j < bases[mp].ninterface; j++) {
	      curry = bases[mp].y[j];
	      if(yboxl/2.0-curry<poswall)
		poswall=yboxl/2.0-curry;//shortest distance in y from the +wall
	      if(curry+yboxl/2.0<negwall)
		negwall=curry+yboxl/2.0;//shortest distance in y from the -wall. >0 inside the box, <0 outside the box.
	     
	      ychg = bases[mp].y[j] - yboxl / 2.0;
	      if (ychg > 0) {
		flagpos = 1;
		if (-ychg < ytot)
		  ytot = -ychg;
	      }
	      if (ychg < -yboxl) {
		flagneg = 1;
		if (-(ychg + yboxl) > ytot)
		  ytot = -(ychg + yboxl);
	      }
	    }
	  }//finished collecting flagpos and flagneg for complex 1
	  for (i = 0; i < s2; i++) {
	    mp = ind_com[k2].plist[i];
	    /*measure each interface to y plane*/
	    for (j = 0; j < bases[mp].ninterface; j++) {
	      curry = bases[mp].y[j];
	      if(yboxl/2.0-curry<poswall)
		poswall=yboxl/2.0-curry;//shortest distance in y from the +wall
	      if(curry+yboxl/2.0<negwall)
		negwall=curry+yboxl/2.0;//shortest distance in y from the -wall. >0 inside the box, <0 outside the box.
	     
	      ychg = bases[mp].y[j] - yboxl / 2.0;
	      if (ychg > 0) {
		flagpos = 1;
		if (-ychg < ytot2)
		  ytot2 = -ychg;
	      }
	      if (ychg < -yboxl) {
		flagneg = 1;
		if (-(ychg + yboxl) > ytot2)
		  ytot2 = -(ychg + yboxl);
	      }
	    }
	  }//finished collecting flagpos and flagneg for complex 2
	  //cout <<"flagneg and flagpos:" <<flagneg<<' '<<flagpos<<endl;
	  //cout <<"poswal and negwall: "<<poswall<<' '<<negwall<<endl;
	  
	  if (flagneg >0 && flagpos>0) {
	    cancel=-1;//current pos spans both
	    cout <<"STICKS OUT BOTH SIDES in Y, CANCEL ASSOCIATION "<<endl;
	    //write_crds_complex(k, ind_com, bases);
	    //write_crds_complex(k2, ind_com, bases);
	  }else if(flagneg>0){
	    //ytot is positive.
	    if(poswall<ytot || poswall<ytot2){
	      cancel=-1;//current pos spans both
	      cout <<"WILL STICK OUT in Y, CANCEL ASSOCIATION "<<endl;
	      // write_crds_complex(k, ind_com, bases);
	      //write_crds_complex(k2, ind_com, bases);
	    }else{
	      //put back in the box. put at edge, rather than bouncing off.
	      if(ytot2>ytot)ytot=ytot2;
	      ind_com[k].ycom +=  ytot;
	      for (i = 0; i < s1; i++) {
		mp = ind_com[k].plist[i];
		bases[mp].ycom +=  ytot;
		  //update interface coords
		  for (j = 0; j < bases[mp].ninterface; j++) {
		    bases[mp].y[j] += ytot;
		  }
	      }
	      ind_com[k2].ycom += ytot;
	      for (i = 0; i < s2; i++) {
		  mp = ind_com[k2].plist[i];
		  bases[mp].ycom +=  ytot;
		  //update interface coords
		  for (j = 0; j < bases[mp].ninterface; j++) {
		    bases[mp].y[j] +=  ytot;
		  }
		}

	    }
	  }else if(flagpos>0){
	    //fpos>0, which means ytot is negative.
	    if(negwall<-ytot || negwall<-ytot2){
	      cancel=-1;//current pos spans both
	      cout <<"WILL STICK OUT in Y, CANCEL ASSOCIATION "<<endl;
	      //write_crds_complex(k, ind_com, bases);
	      //write_crds_complex(k2, ind_com, bases);
	    }else{
	      //put back in the box. put at edge, rather than bouncing off.
	      if(ytot2<ytot)ytot=ytot2;
	      ind_com[k].ycom +=  ytot;
	      for (i = 0; i < s1; i++) {
		mp = ind_com[k].plist[i];
		bases[mp].ycom +=  ytot;
		  //update interface coords
		  for (j = 0; j < bases[mp].ninterface; j++) {
		    bases[mp].y[j] += ytot;
		  }
	      }
	      ind_com[k2].ycom += ytot;
	      for (i = 0; i < s2; i++) {
		  mp = ind_com[k2].plist[i];
		  bases[mp].ycom +=  ytot;
		  //update interface coords
		  for (j = 0; j < bases[mp].ninterface; j++) {
		    bases[mp].y[j] +=  ytot;
		  }
		}

	    }
	  }
	  
	}//flagy>0
	
	if (flagx >1) {
	  
	  
	  cout <<"CHECK BOX SPAN in X. Complex 1 Radius, 2 radius "<<rad<<' '<<rad2<<endl;
	  poswall=plist.xboxl;
	  negwall=plist.xboxl;
	  flagpos = 0;
	  flagneg=0;
	  xtot = 0;
	  xtot2=0;
	  for (i = 0; i < s1; i++) {
	    mp = ind_com[k].plist[i];
	    /*measure each interface to x plane*/
	    for (j = 0; j < bases[mp].ninterface; j++) {
	      currx = bases[mp].x[j];
	      if(xboxl/2.0-currx<poswall)
		poswall=xboxl/2.0-currx;//shortest distance in x from the +wall
	      if(currx+xboxl/2.0<negwall)
		negwall=currx+xboxl/2.0;//shortest distance in x from the -wall. >0 inside the box, <0 outside the box.
	     
	      xchg = bases[mp].x[j] - xboxl / 2.0;
	      if (xchg > 0) {
		flagpos = 1;
		if (-xchg < xtot)
		  xtot = -xchg;
	      }
	      if (xchg < -xboxl) {
		flagneg = 1;
		if (-(xchg + xboxl) > xtot)
		  xtot = -(xchg + xboxl);
	      }
	    }
	  }//finished collecting flagpos and flagneg for complex 1
	  for (i = 0; i < s2; i++) {
	    mp = ind_com[k2].plist[i];
	    /*measure each interface to x plane*/
	    for (j = 0; j < bases[mp].ninterface; j++) {
	      currx = bases[mp].x[j];
	      if(xboxl/2.0-currx<poswall)
		poswall=xboxl/2.0-currx;//shortest distance in x from the +wall
	      if(currx+xboxl/2.0<negwall)
		negwall=currx+xboxl/2.0;//shortest distance in x from the -wall. >0 inside the box, <0 outside the box.
	     
	      xchg = bases[mp].x[j] - xboxl / 2.0;
	      if (xchg > 0) {
		flagpos = 1;
		if (-xchg < xtot2)
		  xtot2 = -xchg;
	      }
	      if (xchg < -xboxl) {
		flagneg = 1;
		if (-(xchg + xboxl) > xtot2)
		  xtot2 = -(xchg + xboxl);
	      }
	    }
	  }//finished collecting flagpos and flagneg for complex 2
	  //  cout <<"flagneg and flagpos:" <<flagneg<<' '<<flagpos<<endl;
	  //cout <<"poswal and negwall: "<<poswall<<' '<<negwall<<endl;
	  if (flagneg >0 && flagpos>0) {
	    cancel=-1;//current pos spans both
	    cout <<"STICKS OUT BOTH SIDES in X, CANCEL ASSOCIATION "<<endl;
	    //write_crds_complex(k, ind_com, bases);
	    //write_crds_complex(k2, ind_com, bases);
	  }else if(flagneg>0){
	    //xtot is positive.
	    if(poswall<xtot || poswall<xtot2){
	      cancel=-1;//current pos spans both
	      cout <<"WILL STICK OUT in X, CANCEL ASSOCIATION "<<endl;
	      //write_crds_complex(k, ind_com, bases);
	      //write_crds_complex(k2, ind_com, bases);
	    }else{
	      //put back in the box. put at edge, rather than bouncing off.
	      if(xtot2>xtot)xtot=xtot2;
	      ind_com[k].xcom +=  xtot;
	      for (i = 0; i < s1; i++) {
		mp = ind_com[k].plist[i];
		bases[mp].xcom +=  xtot;
		  //update interface coords
		  for (j = 0; j < bases[mp].ninterface; j++) {
		    bases[mp].x[j] += xtot;
		  }
	      }
	      ind_com[k2].xcom += xtot;
	      for (i = 0; i < s2; i++) {
		  mp = ind_com[k2].plist[i];
		  bases[mp].xcom +=  xtot;
		  //update interface coords
		  for (j = 0; j < bases[mp].ninterface; j++) {
		    bases[mp].x[j] +=  xtot;
		  }
		}

	    }
	  }else if(flagpos>0){
	    //fpos>0, which means xtot is negative.
	    if(negwall<-xtot || negwall<-xtot2){
	      cancel=-1;//current pos spans both
	      cout <<"WILL STICK OUT in X, CANCEL ASSOCIATION "<<endl;
	      //write_crds_complex(k, ind_com, bases);
	      //write_crds_complex(k2, ind_com, bases);
	    }else{
	      //put back in the box. put at edge, rather than bouncing off.
	      if(xtot2<xtot)xtot=xtot2;
	      ind_com[k].xcom +=  xtot;
	      for (i = 0; i < s1; i++) {
		mp = ind_com[k].plist[i];
		bases[mp].xcom +=  xtot;
		  //update interface coords
		  for (j = 0; j < bases[mp].ninterface; j++) {
		    bases[mp].x[j] += xtot;
		  }
	      }
	      ind_com[k2].xcom += xtot;
	      for (i = 0; i < s2; i++) {
		  mp = ind_com[k2].plist[i];
		  bases[mp].xcom +=  xtot;
		  //update interface coords
		  for (j = 0; j < bases[mp].ninterface; j++) {
		    bases[mp].x[j] +=  xtot;
		  }
		}

	    }
	  }
	}
	if(cancel==-1)
	  cout <<" CANCEL, NEW COMPLEX IS TOO LARGE FOR THE BOX SIZE! DO NOT PERFORM ASSOCIATION. "<<endl;
	return cancel;
}
