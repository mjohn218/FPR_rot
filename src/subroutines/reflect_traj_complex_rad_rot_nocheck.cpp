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

void reflect_traj_complex_rad_rot_nocheck(int p1, Fullmol *bases, Complex *ind_com, double xboxl, double yboxl, double zboxl, double *traj, double *M) {
	int i, j;
	int k = bases[p1].mycomplex;
	int s1 = ind_com[k].mysize;
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
	xchg = currx - xboxl / 2.0;

	if ((xchg + rad) > 0) {

		flagx = 1;
	} else if ((xchg - rad) < -xboxl) {

		flagx += 1;
	}
	ychg = curry - yboxl / 2.0;
	if ((ychg + rad) > 0) {

		flagy = 1;
	} else if ((ychg - rad) < -yboxl) {

		flagy += 1;
	}
	zchg = currz - zboxl / 2.0;
	if ((zchg + rad) > 0) {

		flagz = 1;
	} else if ((zchg - rad) < -zboxl) {

		flagz += 1;
	}

	int mp;

	/*Z is separate to allow the interfaces to approach to the membrane
	 but don't need to test if the entire complex is far enough
	 away from the boundary.
	 */
	
	double dx, dy, dz;
	double x0, y0, z0;
	double dzrot, dxrot, dyrot;
	double row[3];
	int flag;
	if (flagx >0) {
		flag = 0;
		xtot = 0;
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
			
			/*Measure each protein COM to z plane*/
			dx = bases[mp].xcom - x0;
			dy = bases[mp].ycom - y0;
			dz = bases[mp].zcom - z0;
			/*only need z component*/
			//rotate(v, M, v2);
			dxrot = row[0] * dx + row[1] * dy + row[2] * dz;
			currx = x0 + traj[0] + dxrot;
			
			xchg = currx - xboxl / 2.0;
			if (xchg > 0) {
			  flag = 1;
			  if (-xchg < xtot)
			    xtot = -xchg;
			  
			} else if (xchg < -xboxl) {
			  flag = 1;
			  if (-(xchg + xboxl) > xtot)
			    xtot = -(xchg + xboxl);
			}
			/*measure each interface to z plane*/
			for (j = 0; j < bases[mp].ninterface; j++) {
				dx = bases[mp].x[j] - x0;
				dy = bases[mp].y[j] - y0;
				dz = bases[mp].z[j] - z0;
				/*only need z component*/
				//rotate(v, M, v2);
				dxrot = row[0] * dx + row[1] * dy + row[2] * dz;
				currx = x0 + traj[0] + dxrot;

				xchg = currx - xboxl / 2.0;
				if (xchg > 0) {
					flag = 1;
					if (-xchg < xtot)
						xtot = -xchg;

				} else if (xchg < -xboxl) {
					flag = 1;
					if (-(xchg + xboxl) > xtot)
						xtot = -(xchg + xboxl);
				}
			}
		}
		if (flag == 1) {

		  /*Put back inside the box, extended out*/
		  traj[0] += 2.0 * xtot;
		  
		}


	}
	if (flagy >0) {
		flag = 0;
		ytot = 0;
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
			
			ychg = curry - yboxl / 2.0;
			if (ychg > 0) {
			  flag = 1;
			  if (-ychg < ytot)
			    ytot = -ychg;
			  
			} else if (ychg < -yboxl) {
			  flag = 1;
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

				ychg = curry - yboxl / 2.0;
				if (ychg > 0) {
					flag = 1;
					if (-ychg < ytot)
						ytot = -ychg;

				} else if (ychg < -yboxl) {
					flag = 1;
					if (-(ychg + yboxl) > ytot)
						ytot = -(ychg + yboxl);
				}
			}
		}
		if (flag == 1) {

		  /*Put back inside the box*/
		  traj[1] += 2.0 * ytot;
		}
		
	}
	if (flagz >0) {
		flag = 0;
		ztot = 0;
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
			
			/*measure each protein com to z plane*/
			dx = bases[mp].xcom - x0;
			dy = bases[mp].ycom - y0;
			dz = bases[mp].zcom - z0;
			/*only need z component*/
			//rotate(v, M, v2);
			dzrot = row[0] * dx + row[1] * dy + row[2] * dz;
			currz = z0 + traj[2] + dzrot;
			
			zchg = currz - zboxl / 2.0;
			if (zchg > 0) {
			  flag = 1;
			  if (-zchg < ztot)
			    ztot = -zchg;
			  
			} else if (zchg < -zboxl) {
			  flag = 1;
			  if (-(zchg + zboxl) > ztot)
			    ztot = -(zchg + zboxl);
			}
			/*measure each interface to z plane*/
			for (j = 0; j < bases[mp].ninterface; j++) {
			  dx = bases[mp].x[j] - x0;
				dy = bases[mp].y[j] - y0;
				dz = bases[mp].z[j] - z0;
				/*only need z component*/
				//rotate(v, M, v2);
				dzrot = row[0] * dx + row[1] * dy + row[2] * dz;
				currz = z0 + traj[2] + dzrot;

				zchg = currz - zboxl / 2.0;
				if (zchg > 0) {
					flag = 1;
					if (-zchg < ztot)
						ztot = -zchg;

				} else if (zchg < -zboxl) {
					flag = 1;
					if (-(zchg + zboxl) > ztot)
						ztot = -(zchg + zboxl);
				}
			}
		}
		if (flag == 1) {

			/*Put back inside the box*/
			traj[2] += 2.0 * ztot;
		
		}
		
	}

}
