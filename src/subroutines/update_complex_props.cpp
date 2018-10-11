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

void update_complex_props(int c1, Parms plist, Complex *ind_com, Fullmol *bases)
{

  update_one_com_only(c1, ind_com, bases);
  update_radius(c1, ind_com, bases);
  update_diffusion(c1, ind_com, bases);
  update_rot_diffusion(c1, ind_com, bases);
  

}
