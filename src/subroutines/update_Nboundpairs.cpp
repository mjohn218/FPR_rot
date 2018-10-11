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

void update_Nboundpairs(int ptype1, int ptype2, int chg,  Parms plist, int *nBoundPairs)
{
  /*After a reaction occurs, change the count of protein-protein pair bonds */
  /*To avoid storing the same pair in two bins, ptype1 must be <=ptype2*/
  
  int Npro=plist.Nprotypes;
  
  int ind=ptype1*Npro+ptype2;
  if(ptype1>ptype2)ind=ptype2*Npro+ptype1;
  //cout <<"Change Nbound pairs from : "<<nBoundPairs[ind]<<" by: "<<chg<<"  index: "<<ind<<endl;
  nBoundPairs[ind]=nBoundPairs[ind]+chg;
  
}
