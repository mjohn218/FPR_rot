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

void remove_one_prob_all(int p1, int *ncross, int **crosspart, double **prob, int **cross_rxn, int **cross_int) {

	/*If a protein already tried to react and probability was too low,
	 make sure not to try to react it again when it is the partner's turn to
	 test for reactions by setting the reaction probability to zero
	 */

	int pskip;
	int i, j;

	for (i = 0; i < ncross[p1]; i++) {

		pskip = crosspart[p1][i];
		/*for each partner, set it's reaction prob with p1 to zero*/
		for (j = 0; j < ncross[pskip]; j++) {
			if (crosspart[pskip][j] == p1) {
				prob[pskip][j] = 0; //set to zero, that way you won't try to react, but you will diffuse and avoid.

			}
		}

	}

}
