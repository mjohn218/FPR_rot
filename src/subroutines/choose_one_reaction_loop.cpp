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

int choose_one_reaction_loop(int p1, int *ncross, int **crosspart, double **prob, int &ci1, int &ci2, int **cross_rxn, int **cross_int, double irandmax) {

	/*
	  Compare probability of reacting with all partners to random num, if 
	 reaction probability is larger than rnum, select the individual reaction
	 who's probability spans the rnum.
	 
	 In this version, loop over each partner separately, instead of adding probabilities.
	*/
	int p2, pskip;
	int i, j;

	double pmatch;
	double pnow = 0;
	double pmax = 0;
	int flag = 0;
	double rnum;
	double rnum2;
	for (i = 0; i < ncross[p1]; i++) {
	  
	  rnum=1.0*rand_gsl();
	  
	  if (rnum > prob[p1][i]) {
	    flag = 0;
	    /*No reaction occurs, save current partners so can test not overlap*/

	  } else {
	    rnum2 = rnum + rand_gsl() * irandmax; //refine the random number so we don't get exactly zero!
	    if (rnum2 > prob[p1][i]) {
	      flag = 0;
	    } else {
	      /*We will have an association reaction*/
	      flag = 1; //we chose an association reaction
	      ci1 = i;
	      i=ncross[p1]; //Break out of loop, only choose one reaction per protein
	      
	      p2 = crosspart[p1][ci1];
	      pmatch = prob[p1][ci1];
	      for (j = 0; j < ncross[p2]; j++) {
		if (prob[p2][j] == pmatch) {
		  /*it is possible that you had same prob for multiple partners, make sure rxn and partners match*/
		  ci2 = j;
		  if (crosspart[p2][ci2] == p1){
		    if(cross_rxn[p1][ci1]==cross_rxn[p2][ci2])
		      j = ncross[p2]; //break from loop, otherwise the next prob 0 point could be your reactant
		  }
		}
	      }
	      
	    }//will perform association
	  }
	}
	return flag;
}
