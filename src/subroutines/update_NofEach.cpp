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

void update_NofEach(int c1, Parms plist, Complex *ind_com, Fullmol *bases)
{
  int p, mp, ptype;
  for(p=0;p<plist.Nprotypes;p++)
    ind_com[c1].NofEach[p]=0;//reset to zero, after dissociation.
  
  for(p=0;p<ind_com[c1].mysize;p++){
    mp=ind_com[c1].plist[p];
    ptype=bases[mp].protype;
    ind_com[c1].NofEach[ptype]+=1;
  }
  
}
