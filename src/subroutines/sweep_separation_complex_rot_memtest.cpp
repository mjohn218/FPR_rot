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



void sweep_separation_complex_rot_memtest(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, double **trajR, double *M, int **Rlist, int iter) {
	/*
	  In this version, complex k1 is on the membrane. 
	  If both proteins are on the membrane (Dz==0), evaluate only xy displacement, not z.
		  
	  In this _complex_ version, it tests overlap not just for each protein, but for each complex, so all the proteins in a complex, before performing
	  position updates.
	  
	  NEW 2018: IN THIS VERSION, IT DOES NOT ATTEMPT TO SOLVE OVERLAP FOR PROTEINS WITHIN THE SAME COMPLEX, SINCE THEY CANNOT DIFFUSE RELATIVE TO ONE ANOTHER!

	  
	  This routine checks whether protein p1 is overlapping any partners in its reaction
	  zone at its new position that is given by its current position +traj. If it does
	 overlap, the displacement by traj is rejected and a new position for itself and any overlapping partners are selected. Once it
	 no longer overlaps anyone, this protein and its complex are moved and the partners
	 retain their stored new displacements. If a protein has already updated its position (done in sequential order)
	 then it cannot resample a new position, the current protein must still continue to avoid overlapping, however.
	
	 

	*/
  
	int i, j;
	double *M2 = new double[9];
	double *v = new double[3];
	double *v2 = new double[3];
	
	int k1, k2;
	double dx1, dx2;
	double dy1, dy2;
	double dz1, dz2;
	int p2;
	int i1, i2;
	int rxn;
	int iind1, iind2;
	double dr2;
	k1 = bases[p1].mycomplex;
	int pstart = p1;

	int csize = ind_com[k1].mysize;
        int maxrow = 1;
	int c = 0, tp = 0;
  
	for( c = 0; c < csize; c++ )
	{
		tp = ind_com[k1].plist[c];

		if( ncross[tp] > maxrow )
			maxrow = ncross[tp];
	}

//	int row = ncross[p1] + 5; //+5 is buffer


	int ilist[maxrow * csize];
	int olaplist[maxrow * csize];
	int memtest[maxrow *csize];
	/*The sampled displacement for p1 is stored in traj. the position from the
	 previous step is still stored in bases[p1].xcom, etc, and will be updated
	 at the end of this routine*/

	double df1, df2, df3;
	/*figure out i2*/
	
	for (c = 0; c < csize; c++) {
		p1 = ind_com[k1].plist[c];
		for (i = 0; i < ncross[p1]; i++) {
			p2 = cross_part[p1][i];
			k2 = bases[p2].mycomplex;
			if(ind_com[k2].Dz==0){
			  memtest[maxrow*c+i]=1;
			  //cout <<"2D distance in sweep: "<<p1<<" p2: "<<p2<<" complex: "<<k1<<" k2: "<<k2<<" size k1: "<<ind_com[k1].mysize<<" size k2: "<<ind_com[k2].mysize<<" movestat[p1]: "<<movestat[p1]<<" movestat[p2] "<<movestat[p2]<<endl;
			}else
			  memtest[maxrow*c+i]=0;
			i1 = cross_int[p1][i];
			rxn = cross_rxn[p1][i];
			if (Rlist[rxn][0] == i1)
			  i2 = Rlist[rxn][1]; //this is the interface we're looking for on the other protein
			else
			  i2 = Rlist[rxn][0];
			
			ilist[maxrow * c + i] = i2;
		}
	}
	int it = 0;
	int maxit = 50;
	int saveit = 0;
	int t = 0;
	int flag = 0;
	int tsave = 0;
	int flags = 0;
	//   if(csize>1){
	//     for(c=0;c<csize;c++){
	//       p1=ind_com[k1].plist[c];
	//       rxn=cross_rxn[p1][0];
	//       cout <<"complex of 2: "<<p1<<" bound to: "<<bases[p1].partner[0]<<" ncross: "<<ncross[p1]<<" binrdad: "<<bindrad[rxn]<<" first partner: "<<cross_part[p1][0]<<endl;
	//     }
	//   }
	double x0 = ind_com[k1].xcom;
	double y0 = ind_com[k1].ycom;
	double z0 = ind_com[k1].zcom;
	rotationEuler(trajR[k1][0], trajR[k1][1], trajR[k1][2], M); //M is always for protein p1
	reflect_traj_complex_rad_rot(p1, bases, ind_com, plist, traj[k1], trajR[k1], M);
	double x02;
	double y02;
	double z02;


	while (it < maxit) {
		t = 0;
		flag = 0;
		for (c = 0; c < csize; c++) {
			p1 = ind_com[k1].plist[c];
			for (i = 0; i < ncross[p1]; i++) {
				p2 = cross_part[p1][i];
				k2 = bases[p2].mycomplex;
				/*Do not sweep for overlap if proteins are in the same complex, they cannot move relative to one another!
				 */
				if(k1!=k2){
				  i1 = cross_int[p1][i];
				  rxn = cross_rxn[p1][i];
				  i2 = ilist[maxrow * c + i];
				  iind1 = ihome[i1];
				  iind2 = ihome[i2];
				  
				  v[0] = bases[p1].x[iind1] - x0;
				  v[1] = bases[p1].y[iind1] - y0;
				  v[2] = bases[p1].z[iind1] - z0;
				  
				  
				  rotate(v, M, v2);
				  dx1 = x0 + v2[0] + traj[k1][0]; 
				  dy1 = y0 + v2[1] + traj[k1][1]; 
				  dz1 = z0 + v2[2] + traj[k1][2]; 
				  // if(k1==1138 && iter>127000){
				  //   cout <<" In sweep, complex 1138. iter "<<iter << " loop it: "<<it<<endl;
				  //   cout <<"traj:" <<traj[k1][0]<<' '<<traj[k1][1]<<' '<<traj[k1][2]<<endl;
				  //   cout <<"traj:" <<trajR[k1][0]<<' '<<trajR[k1][1]<<' '<<trajR[k1][2]<<endl;
				  // }
				  
				  /*Now complex 2*/
				  
				  rotationEuler(trajR[k2][0], trajR[k2][1], trajR[k2][2], M2);
				  reflect_traj_complex_rad_rot(p2, bases, ind_com, plist, traj[k2], trajR[k2], M2);
				  x02 = ind_com[k2].xcom;
				  y02 = ind_com[k2].ycom;
				  z02 = ind_com[k2].zcom;
				  
				  v[0] = bases[p2].x[iind2] - x02;
				  v[1] = bases[p2].y[iind2] - y02;
				  v[2] = bases[p2].z[iind2] - z02;
				  
				  
				  rotate(v, M2, v2);
				  
				  dx2 = x02 + v2[0] + traj[k2][0];
				  dy2 = y02 + v2[1] + traj[k2][1];
				  dz2 = z02 + v2[2] + traj[k2][2];
				  // if(k2==1138 && iter>127000){
				  //   cout <<" In sweep, complex (k2) 1138. iter "<<iter << " loop it: "<<it<<endl;
				  //   cout <<"traj:" <<traj[k2][0]<<' '<<traj[k2][1]<<' '<<traj[k2][2]<<endl;
				  //   cout <<"traj:" <<trajR[k2][0]<<' '<<trajR[k2][1]<<' '<<trajR[k2][2]<<endl;
				  // }
				  
				  /*separation*/
				  df1 = dx1 - dx2;
				  df2 = dy1 - dy2;
				  df3 = dz1 - dz2;
				  
				  if(memtest[maxrow * c + i]==1){
				    dr2 = df1 * df1 + df2 * df2; //xy displacement only
				    //cout <<" distance in memtest sweep between: "<<p1<<" p2: "<<p2<<" dr2: "<<dr2<<" 3D: "<< df1 * df1 + df2 * df2 + df3 * df3<<endl;
				  }else
				    dr2 = df1 * df1 + df2 * df2 + df3 * df3;
				  
				  if (dr2 < bindrad[rxn] * bindrad[rxn]) {
				    /*reselect positions for protein p2*/
				    // cout<<"reselect: "<<it<<' '<<sqrt(dr2)<<" p1: "<<p1<<" p2: "<<p2<<'\t'<<" dr2: "<<dr2<<endl;
				    olaplist[t] = p2;
				    t++;
				    flag = 1;
				  }
				}//ignore proteins within the same complex.
			}
		}
		/*Now resample positions of p1 and olaplist, if t>0, otherwise no overlap, so
		 break from loop*/
		if (flag == 1) {
			it++;
			traj[k1][0] = sqrt(2.0 * deltat * ind_com[k1].Dx) * GaussV();
			traj[k1][1] = sqrt(2.0 * deltat * ind_com[k1].Dy) * GaussV();
			traj[k1][2] = sqrt(2.0 * deltat * ind_com[k1].Dz) * GaussV();
			trajR[k1][0] = sqrt(2.0 * deltat * ind_com[k1].Drx) * GaussV();
			trajR[k1][1] = sqrt(2.0 * deltat * ind_com[k1].Dry) * GaussV();
			trajR[k1][2] = sqrt(2.0 * deltat * ind_com[k1].Drz) * GaussV();
			
			rotationEuler(trajR[k1][0], trajR[k1][1], trajR[k1][2], M);
			
			reflect_traj_complex_rad_rot(p1, bases, ind_com, plist, traj[k1], trajR[k1], M);
			x0 = ind_com[k1].xcom;
			y0 = ind_com[k1].ycom;
			z0 = ind_com[k1].zcom;
			
			//	  cout<<'p1'<<' '<<traj[k1][0]<<' '<<traj[k1][1]<<' '<<ind_com[k1].Dx<<' '<<ind_com[k1].Dy<<endl;
			for (j = 0; j < t; j++) {
				p2 = olaplist[j];
				k2 = bases[p2].mycomplex;
				if (p2 > pstart && movestat[p2] != 2) {
					/*
					 We loop over proteins sequentially, so earlier proteins have already moved and avoided
					 their neighbors and should not be moved again.
					 These new positions selected for proteins not yet moved will be stored and
					 then used when they test for overlap themselves.
					 */

					/*If p2 just dissociated, also don't try to move again*/
					traj[k2][0] = sqrt(2.0 * deltat * ind_com[k2].Dx) * GaussV();
					traj[k2][1] = sqrt(2.0 * deltat * ind_com[k2].Dy) * GaussV();
					traj[k2][2] = sqrt(2.0 * deltat * ind_com[k2].Dz) * GaussV();
					//	  cout<<'p2'<<' '<<traj[k2][0]<<' '<<traj[k2][1]<<' '<<ind_com[k2].Dx<<' '<<ind_com[k2].Dy<<endl;
					trajR[k2][0] = sqrt(2.0 * deltat * ind_com[k2].Drx) * GaussV();
					trajR[k2][1] = sqrt(2.0 * deltat * ind_com[k2].Dry) * GaussV();
					trajR[k2][2] = sqrt(2.0 * deltat * ind_com[k2].Drz) * GaussV();
					rotationEuler(trajR[k2][0], trajR[k2][1], trajR[k2][2], M2);
					reflect_traj_complex_rad_rot(p2, bases, ind_com, plist, traj[k2], trajR[k2], M2);
				}
			}
			tsave = t;

		} else {
			saveit = it;
			it = maxit; //break from loop
		}

		if (it == maxit - 1) {
		  if (k1 != k2) {
		    cout <<" WARNING ***************************************************** "<<endl;
		    cout << "can't solve overlap: " << p1 << " max cross: " << ncross[p1] << " n overlap: " << tsave << " pro1: " << olaplist[0] << " D: " << ind_com[k1].Dx << " " << ind_com[k2].Dx << " Last separation: "<<sqrt(dr2)<<endl;
		  
		    //write_crds_complex(k1, ind_com, bases);
		    //write_crds_complex(k2, ind_com, bases);
		    for(c=0;c<csize;c++){
		      p1=ind_com[k1].plist[c];
		      //cout <<"complex size: "<<csize<<" cvalue: "<<c<<endl;
		      cout <<"p1: "<<p1<<' '<<" nfree and com; "<<bases[p1].nfree<<' '<<bases[p1].xcom<<' '<<bases[p1].ycom<<' '<<bases[p1].zcom<<'\t';
		      cout<<"traj 1: "<<' '<<traj[k1][0]<<' '<<traj[k1][1]<<' '<<traj[k1][2]<<endl;
		      for (i = 0; i < ncross[p1]; i++) {
			p2 = cross_part[p1][i];
			k2 = bases[p2].mycomplex;
			i1 = cross_int[p1][i];
			rxn = cross_rxn[p1][i];
			i2 = ilist[maxrow * c + i];
			iind1 = ihome[i1];
			iind2 = ihome[i2];
			cout <<"cross num: "<<i<<" i1: "<<i1<<" i2: "<<i2<<" pro: "<<p2<<' '<<" nfree; "<<bases[p2].nfree<<' '<<bases[p2].xcom<<' '<<bases[p2].ycom<<' '<<bases[p2].zcom<<'\t';
			cout<<"traj: "<<' '<<traj[k2][0]<<' '<<traj[k2][1]<<' '<<traj[k2][2]<<endl;
			
		      }
		    }
		  }else{
		    cout <<"Reached max iterations in sweep (memtest) without solving overlap, but proteins are within the same complex. Last separation: "<<sqrt(dr2)<<endl;
		  }
		  //exit(1);
		}
		
	} //end maximum iterations
	  //cout <<"solved its: "<<p1<<" ncross[p1]: "<<ncross[p1]<<" nit "<<saveit<<endl;

	/*
	 if(flags==1){
	 //didn't solve, put molecules at contact
	 dx1=traj[k1][0]+bases[p1].x[iind1];
	 dy1=traj[k1][1]+bases[p1].y[iind1];
	 dz1=traj[k1][2]+bases[p1].z[iind1];

	 dx2=traj[k2][0]+bases[p2].x[iind2];
	 dy2=traj[k2][1]+bases[p2].y[iind2];
	 dz2=traj[k2][2]+bases[p2].z[iind2];

	 df1=dx1-dx2;
	 df2=dy1-dy2;
	 df3=dz1-dz2;
	 df1-=plist.xboxl*round(df1/plist.xboxl);
	 df2-=plist.yboxl*round(df2/plist.yboxl);
	 df3-=plist.zboxl*round(df3/plist.zboxl);

	 dr2=df1*df1+df2*df2+df3*df3;
	 double stretch=bindrad[rxn]/sqrt(dr2)+1E-4*rand_gsl();
	 traj[k1][0]+=df1*(stretch-1.0);
	 traj[k1][1]+=df2*(stretch-1.0);
	 traj[k1][2]+=df3*(stretch-1.0);

	 }
	 */

	int s1 = ind_com[k1].mysize;

	dx1 = traj[k1][0];
	dy1 = traj[k1][1];
	dz1 = traj[k1][2];
	int mp;
	ind_com[k1].xcom += dx1;
	ind_com[k1].ycom += dy1;
	ind_com[k1].zcom += dz1;
	
	
	for (i = 0; i < s1; i++) {
		mp = ind_com[k1].plist[i];
		movestat[mp] = 2; //physically updated position of these proteins
		v[0] = bases[mp].xcom - x0;
		v[1] = bases[mp].ycom - y0;
		v[2] = bases[mp].zcom - z0;
		

		rotate(v, M, v2); //M is one used in last step with no overlap.
		bases[mp].xcom = x0 + v2[0] + traj[k1][0];
		bases[mp].ycom = y0 + v2[1] + traj[k1][1];
		bases[mp].zcom = z0 + v2[2] + traj[k1][2];

		
		//update interface coords
		for (j = 0; j < bases[mp].ninterface; j++) {
		  v[0] = bases[mp].x[j] - x0;
		  v[1] = bases[mp].y[j] - y0;
		  v[2] = bases[mp].z[j] - z0;
		  

			rotate(v, M, v2);
			bases[mp].x[j] = x0 + v2[0] + traj[k1][0];
			bases[mp].y[j] = y0 + v2[1] + traj[k1][1];
			bases[mp].z[j] = z0 + v2[2] + traj[k1][2];
			
			
		}
	}
	/*Reset displacements to zero so distance is measured to your current
	  updated position that won't change again this turn
	*/
	
	traj[k1][0] = 0;
	traj[k1][1] = 0;
	traj[k1][2] = 0;
	trajR[k1][0] = 0;
	trajR[k1][1] = 0;
	trajR[k1][2] = 0;
	delete[] M2;
	delete[] v;
	delete[] v2;
	
}
