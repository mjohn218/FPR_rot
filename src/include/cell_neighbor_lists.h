#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>


using namespace std;



void cell_neighbor_listPBCCELL(int Nx, int Ny, int Nz, int maxnbor, int *nbor, int *nborrev);
void cell_neighbor_listPBC(int Nx, int Ny, int Nz, int maxnbor, int *nbor, int *nborrev);
void cell_neighbor_list(int Nx, int Ny, int Nz, int maxnbor, int *nbor, int *Nnbor);
