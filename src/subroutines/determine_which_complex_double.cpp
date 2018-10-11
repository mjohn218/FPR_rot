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

int determine_which_complex_double(int p1, int p2, Complex *ind_com, int c1, int c2, Fullmol *bases, int **myrxn, int **Rlist, int *ihome, int *p_home, Parms &plist) {

  int t1=0;
  int t2=0;
  int k;
  int j, mp;
  int flag;
  int s1=ind_com[c1].mysize;
  int recheck[s1];
  int mycom[s1];
  int rsit[s1];
  int origlist[s1];
  int r=0;
  int cancel=0;

  //upcoming vars
  int nb;
  int itmp;
  int mutmp;
  int i1tmp;
  int iitmp;
  int i2tmp;
  int iitmp2;
  int ppart;
  
  int nre, it;
  int myj;
  int fmatch, cmatch;
  int ptemp, mix, imix;
  int maxit=10000;
  
  /*Here establish which complex all the other proteins are in*/
  //  if(s1>2){
  for(j=0;j<s1;j++){
    mycom[j]=-1;
    rsit[j]=-1;
  }
  maxit=s1*maxit;
  //  cout <<"Size of complex: "<<s1<<" maxit: "<<maxit<<" final complexes: "<<c1<<'\t'<<c2<<endl;

  //otherwise ther's only 2 proteins in this complex
  /*does p1 bind p2 more than once?*/
  nb=bases[p1].nbnd;
  int no=0;
  int dub=0;
  /*THE USE OF THE SINGLE DUB VARIABLE WILL ONLY WORK IF A PROTEIN CAN ONLY HAVE A SINGLE REPEATED BOUND STATE
    OTHERWISE A DUB VARIABLE IS NEEDED FOR EACH BOUND STATE.
   */
  for(k=0;k<nb;k++){
    itmp=bases[p1].bndlist[k];
    mutmp=myrxn[itmp][0];//each bound complex is only in a single reaction
    i1tmp=Rlist[mutmp][1];
    i2tmp=Rlist[mutmp][2];
    /*if both are clathrin, this won't work*/
    //if(p_home[i1tmp]!=bases[p1].protype)
    //  i1tmp=Rlist[mutmp][2];
    iitmp=ihome[i1tmp];
    iitmp2=ihome[i2tmp];
    /*TO GET THE PARTNER PPART OF P1 FOR THIS BOUND COMPLEX*/
    if(bases[p1].istatus[iitmp]==itmp){
      if(iitmp!=iitmp2){
	if(bases[p1].istatus[iitmp2]==itmp){
	  /*if two interfaces have this bndstate, go through both
	    use !=0 because dub can go up to 2 if the middle interface
	    is bound to the same type of interface
	  */
	  if(dub==1)
	    iitmp=iitmp2;
	  dub++;
	}
      }
    }else{
      iitmp=iitmp2;
    }
    ppart=bases[p1].partner[iitmp];
    if(ppart==p2)
      no++;
  }
  if(no>1){
    /*Do not change the sizes of any complexes, just break this bond*/
    cout <<"p1 and p2 doubly bound ! "<<endl;
    bases[p1].mycomplex=c1;
    bases[p2].mycomplex=c1;
    cancel=1;
  }else{
    
    for(j=0;j<s1;j++){
      /*measure the distance between other proteins in the complex*/
      /*test if they are closer to protein i or partner
       *and then move their coordinates accordingly
       */
      mp=ind_com[c1].plist[j];
      origlist[j]=mp;
      //cout <<"Protein in complex: "<<mp<<endl;
      if(mp!=p1 &&mp!=p2){
	
	flag=0;

	/*mp should have all of its binding partners listed*/
	 nb=bases[mp].nbnd;
	 // cout <<"protein: "<<mp<<" nbound: "<<nb<<endl;
	 dub=0;
	 for(k=0;k<nb;k++){
	   itmp=bases[mp].bndlist[k];
	   mutmp=myrxn[itmp][0];
	   i1tmp=Rlist[mutmp][1];
	   i2tmp=Rlist[mutmp][2];
	   //if(p_home[i1tmp]!=bases[mp].protype)
	   //i1tmp=Rlist[mutmp][2];
	   iitmp=ihome[i1tmp];
	   iitmp2=ihome[i2tmp];
	   /*TO GET THE PARTNER PPART OF P1 FOR THIS BOUND COMPLEX*/
	   if(bases[mp].istatus[iitmp]==itmp){
	     if(iitmp!=iitmp2){
	       if(bases[mp].istatus[iitmp2]==itmp){
		 //cout <<"Repeated bound state on : "<<iitmp<<" and " <<iitmp2<<" dub: "<<dub<<" partners: "<<bases[mp].partner[iitmp]<<'\t'<<bases[mp].partner[iitmp2]<<" curr partner: "<<'\t';
		 /*if two interfaces have this bndstate, go through both*/
		 if(dub==1)
		   iitmp=iitmp2;
		 dub++;
		 //cout<<bases[mp].partner[iitmp]<<endl;
	       }
	     }
	   }else{
	     iitmp=iitmp2;
	   }
	   
	  ppart=bases[mp].partner[iitmp];
	  //cout <<"Protein: "<<mp<<" nbnd: "<<nb<<" bndspecie "<<itmp<<" rxn: "<<mutmp<<" iface: "<<i1tmp<<" ihome: "<<iitmp<<" partner: "<<ppart<<endl;
	  if(ppart==p1)
	    flag=1;
	  if(ppart==p2)
	    flag=-1;
	}//end loop k
	
	if(flag==1){
	  //bound with protein p1
	  //cout <<"protein: "<<mp<< " bound to p1: "<<p1<<" nbnd: "<<nb<<'\t';
	  // for(int hh=0;hh<nb;hh++)
	  //cout <<bases[mp].bndlist[hh]<<'\t';
 	  //cout <<endl;
	  ind_com[c1].plist[t1]=mp;//replace [0] with mp so ok, next is current
	  bases[mp].mycomplex=c1;
	  /*move the coordinates of all the interfaces as well*/
	  t1++;
	  mycom[j]=c1;
	}else if(flag==-1){
	  // cout <<"Flag=: "<<flag<<" protein: "<<mp<<" bound to p2: "<<p2<<" nbound: "<<nb<<'\t';
// 	   for(int hh=0;hh<nb;hh++)
//  	    cout <<bases[mp].bndlist[hh]<<'\t';
//  	  cout <<endl;
  
	  //this protein belongs to complex of partner
	  ind_com[c2].plist[t2]=mp;
	  t2++;
	  bases[mp].mycomplex=c2;
	  mycom[j]=c2;
	}else{
	  /*mp is not bound to p1 or to p2, check the other proteins in this complex, and also check whether this
	    complex forms a closed loop*/
	  recheck[r]=mp;
	  rsit[r]=j;
	  r++;
	  
	} //end checking who the other protein in the complex is bound to
	
      }
    }//done looping over all proteins in the original complex
    
    /*now loop over all proteins that need to be rechecked*/
    /*Since two proteins can be bound through multiple interfaces, one also needs to 
      check on unbinding if the proteins that unbound are still in fact bound*/
    
    
    nre=r;
    it=0;
    for(r=0;r<nre;r++){
      mp=recheck[r];
      myj=rsit[r];
      nb=bases[mp].nbnd;
      fmatch=0;
      dub=0;
      it++;
      for(k=0;k<nb;k++){
	itmp=bases[mp].bndlist[k];
	mutmp=myrxn[itmp][0];
	i1tmp=Rlist[mutmp][1];
	i2tmp=Rlist[mutmp][2];
	//if(p_home[i1tmp]!=bases[mp].protype)
	//i1tmp=Rlist[mutmp][2];
	iitmp=ihome[i1tmp];
	iitmp2=ihome[i2tmp];
	/*TO GET THE PARTNER PPART OF P1 FOR THIS BOUND COMPLEX*/
	if(bases[mp].istatus[iitmp]==itmp){
	  if(iitmp!=iitmp2){
	    if(bases[mp].istatus[iitmp2]==itmp){
	      /*if two interfaces have this bndstate, go through both
		same bound state occurs in self-binding
	      */
	      //cout <<"Repeated bound state on : "<<iitmp<<" and " <<iitmp2<<" dub: "<<dub<<" partners: "<<bases[mp].partner[iitmp]<<'\t'<<bases[mp].partner[iitmp2]<<" curr partner: "<<'\t';
	      if(dub==1)
		iitmp=iitmp2;
	      
	      dub++;
	      //cout<<bases[mp].partner[iitmp]<<endl;
	    }
	  }
	}else{
	  iitmp=iitmp2;
	}
	
	ppart=bases[mp].partner[iitmp];
	//cout <<"Protein rechecking: "<<mp<<" nbnd: "<<nb<<" bndspecie "<<itmp<<" rxn: "<<mutmp<<" iface: "<<i1tmp<<" ihome: "<<iitmp<<" partner: "<<ppart<<endl;
	
	/*now loop over hte proteins in this complex to see who you are bound to */
	for(j=0;j<s1;j++){
	  ptemp=origlist[j];
	  if(ptemp!=mp){
	    /*you don't bind to yourself*/
	    if(ppart==ptemp){
		/*either you know this proteins complex or not*/
	      if(mycom[j]==c1){
		ind_com[c1].plist[t1]=mp;//replace [0] with mp so ok, next is current
		bases[mp].mycomplex=c1;
		/*move the coordinates of all the interfaces as well*/
		t1++;
		mycom[myj]=c1;
		//cout <<"rechecked: protein: "<<mp<<" bound to partner: "<<ppart<<" complex c1: "<<c1<<endl;
		//break out
		j=s1;
		k=nb;
		fmatch=1;
	      }else if(mycom[j]==c2){
		ind_com[c2].plist[t2]=mp;//replace [0] with mp so ok, next is current
		bases[mp].mycomplex=c2;
		/*move the coordinates of all the interfaces as well*/
		t2++;
		mycom[myj]=c2;
		//cout <<"rechecked: protein: "<<mp<<" bound to partner: "<<ppart<<" complex c2: "<<c2<<endl;
		//break out
		j=s1;
		k=nb;
		fmatch=1;
	      }
	    }
	  }//ignore yourself
	}//done checking if your partner is part of the complex
	
	
      }//end looping over mp's partners
      if(it%100000==0)
	cout <<"fmatch : "<<fmatch<<" it: "<<it<<" nrecheck: "<<nre<<" curr r: "<<r<<endl;
      if(fmatch==0){
	/*need to keep checking! put this protein last in line*/
	mix=recheck[r];
	imix=rsit[r];
	/*reorder all the ones that need rechecking*/
	for(int rr=r;rr<nre-1;rr++){
	  recheck[rr]=recheck[rr+1];
	  rsit[rr]=rsit[rr+1];
	}
	recheck[nre-1]=mix;
	rsit[nre-1]=imix;
	r--;
	
      }
      if(it>maxit){
	cerr <<"stuck in infinite loop"  <<endl;
	cout <<"stuck in infinite loop"  <<endl;
	for(j=0;j<s1;j++){
	  ptemp=origlist[j];
	  cout <<"original protein: "<<ptemp<<" mycomplex:" <<mycom[j]<<" "<<bases[ptemp].mycomplex<<endl;
	}
	int rs=r+1;
	for(r=rs;r<nre-1;r++){
	  mp=recheck[r];
	  myj=rsit[r];
	  nb=bases[mp].nbnd;
	  fmatch=0;
	  dub=0;
	  //it++;
	  for(k=0;k<nb;k++){
	    itmp=bases[mp].bndlist[k];
	    mutmp=myrxn[itmp][0];
	    i1tmp=Rlist[mutmp][1];
	    i2tmp=Rlist[mutmp][2];
	    iitmp=ihome[i1tmp];
	    iitmp2=ihome[i2tmp];
	    /*TO GET THE PARTNER PPART OF P1 FOR THIS BOUND COMPLEX*/
	    if(bases[mp].istatus[iitmp]==itmp){
	      if(iitmp!=iitmp2){
		if(bases[mp].istatus[iitmp2]==itmp){
		  /*if two interfaces have this bndstate, go through both*/
		  cout <<"Repeated bound state on : "<<iitmp<<" and " <<iitmp2<<" dub: "<<dub<<" partners: "<<bases[mp].partner[iitmp]<<'\t'<<bases[mp].partner[iitmp2]<<" curr partner: "<<'\t';
		  if(dub==1)
		    iitmp=iitmp2;
		
		  dub++;
		}
	      }
	    }else{
	      iitmp=iitmp2;
	    }
	    
	    ppart=bases[mp].partner[iitmp];
	    cout <<"Final dump, Protein rechecking: "<<mp<<" nbnd: "<<nb<<" bndspecie "<<itmp<<" rxn: "<<mutmp<<" iface: "<<i1tmp<<" ihome: "<<iitmp<<" partner: "<<ppart<<endl;
	    
	    /*now loop over hte proteins in this complex to see who you are bound to */
	    for(j=0;j<s1;j++){
	      ptemp=origlist[j];
	      if(ptemp!=mp){
		/*you don't bind to yourself*/
		if(ppart==ptemp){
		  /*either you know this proteins complex or not*/
		  if(mycom[j]==c1){
		    ind_com[c1].plist[t1]=mp;//replace [0] with mp so ok, next is current
		    bases[mp].mycomplex=c1;
		    /*move the coordinates of all the interfaces as well*/
		    t1++;
		    mycom[myj]=c1;
		    cout <<"rechecked: protein: "<<mp<<" bound to partner: "<<ppart<<" complex c1: "<<c1<<endl;
		    //break out
		    j=s1;
		    k=nb;
		    fmatch=1;
		  }else if(mycom[j]==c2){
		    ind_com[c2].plist[t2]=mp;//replace [0] with mp so ok, next is current
		    bases[mp].mycomplex=c2;
		    /*move the coordinates of all the interfaces as well*/
		    t2++;
		    mycom[myj]=c2;
		    cout <<"rechecked: protein: "<<mp<<" bound to partner: "<<ppart<<" complex c2: "<<c2<<endl;
		    //break out
		    j=s1;
		    k=nb;
		    fmatch=1;
		  }else
		    cout<<"evaluated protein: "<<mp<<" to partner : "<<ppart<<" but no complex yet for either "<<endl;
		}
	      }//ignore yourself
	    }//done checking if your partner is part of the complex
	    
	    
	  }//end looping over mp's partners
	}  
	exit(1);
      }
      
    }//end establishing other binding partners
    //cout <<"finshed rechecking "<<endl;
    
    /*Check for closed loops by seeing if all the proteins you are bound to are still in your same complex*/
    for(j=0;j<s1;j++){
      mp=origlist[j];
      if(mp!=p1 &&mp!=p2){
	cmatch=bases[mp].mycomplex;
	/*mp should have all of its binding partners listed*/
	nb=bases[mp].nbnd;
	dub=0;
	for(k=0;k<nb;k++){
	  itmp=bases[mp].bndlist[k];
	  mutmp=myrxn[itmp][0];
	  i1tmp=Rlist[mutmp][1];
	  i2tmp=Rlist[mutmp][2];
	  //if(p_home[i1tmp]!=bases[mp].protype)
	  //i1tmp=Rlist[mutmp][2];
	  iitmp=ihome[i1tmp];
	  iitmp2=ihome[i2tmp];
	  /*TO GET THE PARTNER PPART OF P1 FOR THIS BOUND COMPLEX*/
	  if(bases[mp].istatus[iitmp]==itmp){
	    if(iitmp!=iitmp2){
	      if(bases[mp].istatus[iitmp2]==itmp){
		/*if two interfaces have this bndstate, go through both*/
		if(dub==1)
		  iitmp=iitmp2;
		
		dub++;
	      }
	    }
	  }else{
	    iitmp=iitmp2;
	  }
	  
	  ppart=bases[mp].partner[iitmp];
	  //cout <<"Protein: "<<mp<<" nbnd: "<<nb<<" bndspecie "<<itmp<<" rxn: "<<mutmp<<" iface: "<<i1tmp<<" ihome: "<<iitmp<<" partner: "<<ppart<<endl;
	  if(bases[ppart].mycomplex!=cmatch){
	    cout <<" closed loop, proteins "<<mp<<" and " <<ppart<<" list separate complexes, but still perform dissociation!!  "<<endl;
	    plist.nloop--;
	    cancel=1;
	    k=nb;
	    j=s1;
	  }
	}//end loop k
      }
    }//end looping over all proteins in original complex
    //cout <<"no closed loop "<<endl;
  
    //  }
    ind_com[c1].plist[t1]=p1;
    ind_com[c2].plist[t2]=p2;
    ind_com[c1].mysize=t1+1;//just added in these last proteins, p1 and p2
    ind_com[c2].mysize=t2+1;
    
    bases[p1].mycomplex=c1;
    bases[p2].mycomplex=c2;
  }
  return cancel;

}
