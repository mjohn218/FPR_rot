#include "cell_neighbor_lists.h"

void cell_neighbor_list(int Nx, int Ny, int Nz, int maxnbor, int *nbor, int *Nnbor) {
	int i, j, k;
	int mybin;
	int Ncell = Nx * Ny * Nz;
	int c = 0;
	for (i = 0; i < Ncell; i++)
		Nnbor[i] = 0; //with reflecting boundaries, not all cells will have 13 neighbors.

	for (k = 0; k < Nz; k++) {
		for (j = 0; j < Ny; j++) {
			for (i = 0; i < Nx; i++) {

				/*For each cell figure out its 13 neighbors that are ~forward and up*/
				if (i < (Nx - 1)) {
					//count all the plus x boxes
					mybin = (i + 1) + j * Nx + k * (Nx * Ny);
					nbor[c * maxnbor + Nnbor[c]] = mybin;
					Nnbor[c]++;
					if (j < (Ny - 1)) {
						mybin = (i + 1) + (j + 1) * Nx + k * (Nx * Ny);
						nbor[c * maxnbor + Nnbor[c]] = mybin;
						Nnbor[c]++;
						if (k < (Nz - 1)) {
							mybin = (i + 1) + (j + 1) * Nx + (k + 1) * (Nx * Ny);
							nbor[c * maxnbor + Nnbor[c]] = mybin;
							Nnbor[c]++;

						}

					}
					if (k < (Nz - 1)) {
						mybin = (i + 1) + j * Nx + (k + 1) * (Nx * Ny);
						nbor[c * maxnbor + Nnbor[c]] = mybin;
						Nnbor[c]++;
					}
				}
				if (j < (Ny - 1)) {
					mybin = i + (j + 1) * Nx + k * (Nx * Ny);
					nbor[c * maxnbor + Nnbor[c]] = mybin;
					Nnbor[c]++;
					if (i > 0) {
						mybin = (i - 1) + (j + 1) * Nx + k * (Nx * Ny);
						nbor[c * maxnbor + Nnbor[c]] = mybin;
						Nnbor[c]++;
						if (k < (Nz - 1)) {
							mybin = (i - 1) + (j + 1) * Nx + (k + 1) * (Nx * Ny);
							nbor[c * maxnbor + Nnbor[c]] = mybin;
							Nnbor[c]++;
						}
					}
					if (k < (Nz - 1)) {
						mybin = i + (j + 1) * Nx + (k + 1) * (Nx * Ny);
						nbor[c * maxnbor + Nnbor[c]] = mybin;
						Nnbor[c]++;

					}
				}
				if (k < (Nz - 1)) {
					mybin = i + j * Nx + (k + 1) * (Nx * Ny);
					nbor[c * maxnbor + Nnbor[c]] = mybin;
					Nnbor[c]++;
					if (i > 0) {
						mybin = i - 1 + j * Nx + (k + 1) * (Nx * Ny);
						nbor[c * maxnbor + Nnbor[c]] = mybin;
						Nnbor[c]++;
						if (j > 0) {
							mybin = i - 1 + (j - 1) * Nx + (k + 1) * (Nx * Ny);
							nbor[c * maxnbor + Nnbor[c]] = mybin;
							Nnbor[c]++;
						}
					}
					if (j > 0) {
						mybin = i + (j - 1) * Nx + (k + 1) * (Nx * Ny);
						nbor[c * maxnbor + Nnbor[c]] = mybin;
						Nnbor[c]++;
						if (i < (Nx - 1)) {
							mybin = i + 1 + (j - 1) * Nx + (k + 1) * (Nx * Ny);
							nbor[c * maxnbor + Nnbor[c]] = mybin;
							Nnbor[c]++;
						}

					}

				}
				// cout <<"Cell: "<<c<<" Num neighbors: "<<Nnbor[c]<<'\t';
// 	for(int sn=0;sn<Nnbor[c];sn++)
// 	  cout <<nbor[c*maxnbor+sn]<<" ";
// 	cout <<endl;

				c++;

			} //end looping over z cells
		} //end looping over y cells
	} //end looping over x cells

}
