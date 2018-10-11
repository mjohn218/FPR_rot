/*
 Implement particle-based Reaction diffusion
 algorithm with trajectory reweighting.

 This general version allows any number of different
 particle types and reactions, with rigid body
 motion that includes translation and rotation.

 Periodic boundary conditions are enforced at boundaries.

 */
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>

#include "reactions.h"
#include "rand_gsl.h"
#include "md_timer.h"
#include "GF_calls.h"
#include "utility_calls.h"
#include "vector_rot_calls.h"
#include "cell_neighbor_lists.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_bessel.h>
#include <2Drelated.h>
#include "Faddeeva.hh"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <iomanip>
#include <vector>
#include <string>
#include "evaluate_binding.h"


using namespace std;

struct MD_Timer totaltime;
struct MD_Timer bimoltime;

int main(int argc, char *argv[]) {

	/*Define the reactants, and all the species involved in the reactions
	 *this includes all the products, including the misbinding products
	 */
  int i, idd, j, n, k, recruitmentflag, crosscounter;
	recruitmentflag = 0;
	timeval tim;
	gettimeofday(&tim, 0);
	double t1 = tim.tv_sec + tim.tv_usec;

	vector<gsl_matrix *> contsur;
	vector<gsl_matrix *> contnorm;
	vector<gsl_matrix *> contpir;
	std::vector<std::string> infofilenames;

	int MAXALWTBL = 1000;
	cout << "Maximum number of unique 2D reactions allowed: " << MAXALWTBL << endl;
	double *TBLID=new double[MAXALWTBL*2]; //[type#][0]:ka, [tNx * Ny * (k + 1)ype#][1]:Dtot
	double Dtemp;
	double ktemp;
	int DDtableindex = 0;
	int tableexistflag = 0;
	int uniquetableindex = 0;

	int seed = int(t1);
	//seed=1507126727;
	cout << "seed: " << seed << endl;
	srand_gsl(seed);
	double randmax = pow(2.0, 32);
	double irandmax = 1.0 / randmax;

	
	char fname[100];
	char fnameComplexXYZ[100];
	char fnameProXYZ[100];

	Parms plist;
	plist.restart = 0; //in case you don't read it in.
	plist.pclath = -1; //in case you don't read it in, no clath
	plist.nloop=0;
	
	initialize_timer(&totaltime);
	initialize_timer(&bimoltime);
	start_timer(&totaltime);

	
	/*Input files*/
	ifstream parmfile(argv[1]);
	ifstream rxnfile(argv[2]);
	
	/*Read in parameters*/
	read_parmsN(parmfile, plist,infofilenames); //  write_parms(plist);

	
	int Nprotypes = plist.Nprotypes; //total distinct protein types, so 9
	int Nifaces = plist.Nifaces; //this is the number of interfaces
	int Nrxn = plist.Nrxn;
	int Nspecies = plist.Nspecies; //this will include product species
	int *Ncopy = new int[Nprotypes];

	int Ntotalmol = 0;
	plist.Natom = 0;
	int ntmp;
	
	/*Get copy numbers from Info Files*/
	for (i = 0; i < Nprotypes; i++) {
	  
		ifstream numfile(infofilenames[i].c_str());
		numfile.ignore(400, '\n');
		numfile >> Ncopy[i];
		Ntotalmol += Ncopy[i];
		numfile.close();
	}
	cout << "Ntotal mols: " << Ntotalmol << endl;
	plist.Ntotalmol = Ntotalmol;


	cout <<"Protein type: "<<plist.Nprotypes<<endl;
	for(i=0;i<Nprotypes;i++){
	  cout <<"type: "<<i<<" name length:" <<infofilenames[i].length()<<" 3letname: "<<infofilenames[i].substr(0,3)<<endl;
	}
	/*Declare protein vectors*/
	Fullmol *bases = new Fullmol[Ntotalmol]; //contains information on each protein in the full system
	Complex *ind_com = new Complex[Ntotalmol]; //contains information on each complex
	int *numpartners = new int[Nifaces]; //this should account for all free interfaces
	int **Speclist = new int*[Nifaces];

	for (i = 0; i < Nifaces; i++)
		Speclist[i] = new int[MAXPRTNER];

	Protein *wholep = new Protein[Nprotypes];
	int *p_home = new int[Nifaces]; //this reverses and tells you what protein a given interface belongs to
	int *i_home = new int[Nspecies]; //for both free and bound states, what index are you on the protein's list
	double *bindrad = new double[Nrxn]; //binding or unbinding radius for each reaction
	int *Ncoup = new int[Nrxn]; //list of reactions coupled to this one
	int **mycoupled = new int*[Nrxn];
	for (i = 0; i < Nrxn; i++)
		mycoupled[i] = new int[MAXRXN];
	/*The number of reactions is fixed and all the same reactions are possible in each spatial cell*/
	double *kr = new double[Nrxn]; //reaction rate (with dimensions)

	int **Rlist = new int*[Nrxn]; //the identity of the species in the reaction
	int *Npart = new int[Nrxn]; //The number of participant species in a reaction
	int **Del = new int*[Nrxn]; //The coeffiecients of the participants in the reaction

	int currentnumberofmolectypes = 0;
	int **molectypes; //define and allocate list of unique molecule types
	molectypes = new int*[MAXNUMCOMPLEXTYPES];
	for (i = 0; i < MAXNUMCOMPLEXTYPES; i++)
		molectypes[i] = new int[2+2*Nprotypes];
	for(i=0;i<MAXNUMCOMPLEXTYPES;i++){
	  for(j=0;j<2+2*Nprotypes;j++){
			molectypes[i][j]=0;
		}
	}

	/*Do not assume all reactions possible
	 *so instead, we have to figure out a way to read them in!
	 */
	int maxrctant = 5;
	for (i = 0; i < Nrxn; i++) {
		Rlist[i] = new int[maxrctant];
		Del[i] = new int[maxrctant];
	}
	int *Nmyrxn = new int[Nspecies];
	int **myrxn = new int*[Nspecies];
	for (i = 0; i < Nspecies; i++)
		myrxn[i] = new int[MAXRXN];
	int *cntrxn = new int[3]; //nfree, nbnd, nmut
	int *freelist = new int[Nspecies];
	int *bndlist = new int[Nspecies];
	int *zlist = new int[Nspecies];

	int *rxtype = new int[plist.Nrxn];
        double*** coordcont = new double**[Nprotypes];
        for(i=0; i < Nprotypes; i++){

    	    coordcont[i] = new double*[Nifaces+1];

	        for(j=0; j < Nifaces+1; j++){
	        	coordcont[i][j] = new double[3];

	            for(k=0; k < 3; k++){
	            	coordcont[i][j][k]= 0.0;
	            }
	        }
    	}
	
	int howmanylipids=0;
	
	
	/*Fill in all protein and reaction arrays with input data*/
	read_inputs(rxnfile, bases, Ncopy, p_home, wholep, numpartners, Speclist,  plist, Rlist, Ncoup, mycoupled, Nmyrxn, myrxn, bindrad, kr, cntrxn, freelist, bndlist, zlist, rxtype, i_home, infofilenames, coordcont, howmanylipids);

	/*Generate an initial configuration, prevent overlap of binding partners*/
	generate_initial_crdsNoPBC(plist, bases, Ncopy, ind_com, bindrad, wholep,Rlist,rxtype,p_home,howmanylipids,coordcont);
	
	/*Generate PSF file for reading coordinates into VMD*/
	gen_psf_system(plist, wholep, Ncopy, infofilenames);


	
	int Ntotsite = 0;
	plist.Natomwrite = 0;
	for (i = 0; i < Nprotypes; i++) {
		ntmp = wholep[i].ninterface + 1;
		cout <<"Ninterfaces per protein: "<<i<<" wholep[i].ninterface: "<<wholep[i].ninterface<<endl;
		plist.Natom += Ncopy[i] * ntmp;
		Ntotsite += Ncopy[i] * wholep[i].ninterface;
		plist.Natomwrite += Ncopy[i] * wholep[i].nint_write;
	}
	cout << "N atoms: " << plist.Natom << " Nsites, not COM: " << Ntotsite << " Natoms to write out: " << plist.Natomwrite << endl;
	int t = 0;

	
	/*Change rates to ka and kb, rather than kon and koff
	 unless you are reading in values of ka and kb
	 */
	
	//check_reactions(plist, bases, numpartners, Speclist, Nmyrxn, Rlist, myrxn, Nspecies);

	
	string *names = new string[Nprotypes];
	double *savecrds = new double[Ntotalmol * 3]; //for x, y, z

	copy_crds(Ntotalmol, bases, savecrds);

	/*Print out specific reactions*/
	cout << "Print specific interaction network " << endl;
	int ncomplex = 0;
	for (i = 0; i < Nifaces; i++) {
		cout << i << '\t';
		for (j = 0; j < numpartners[i]; j++) {
			cout << Speclist[i][j] << '\t';
			ncomplex++;
		}
		cout << endl;
	}
	ncomplex /= 2;
	plist.nspec_complex = ncomplex;

	int ind, r1, m;
	int begin, end;

	double tau;

	double rnum;
	double rnum2, rnum3;
	int nfaces = 6; //for a cubic volume
	int direction, neighbor;

	int rxn;

	double curr_time = 0;

	int checkpoint = 10000000; /*How often to write out full current solution*/
	int stepwrite = 100; /*How often to write out species numbers*/
	char fnamemid[100];

	/*G(r) stuff*/
	double delr = 0.1;
	int nbins = 20;
	cout << "nbins: " << nbins << " Ntotal mol: " << Ntotalmol << " max g(r): " << nbins * delr << endl;
	double **gr = new double*[Nprotypes * Nprotypes];
	for (i = 0; i < Nprotypes * Nprotypes; i++)
		gr[i] = new double[nbins];

	int zeroctr=0;
	for(i=0;i<Nrxn;i++){
		if(kr[i]==0){
			zeroctr++;
		}
	}


	char tname[100];
	sprintf(tname, "restartsW.out");

	int newsizeoverlap = int(ceil(MAXOVERLAP*(1+(double)zeroctr/(double)Nrxn)));
	double **probvec = new double*[Ntotalmol];
	int **crosspart = new int*[Ntotalmol]; //Index of the protein partner in this reaction
	int **crossint = new int*[Ntotalmol]; //index of the interface species in this reaction
	int **cross_rxn = new int*[Ntotalmol]; //index of the reaction number of this reaction
	for (i = 0; i < Ntotalmol; i++) {
		probvec[i] = new double[newsizeoverlap];
		crosspart[i] = new int[newsizeoverlap];
		crossint[i] = new int[newsizeoverlap];
		cross_rxn[i] = new int[newsizeoverlap];
	}

	int *ncross = new int[Ntotalmol];
	int *ncrosscom=new int[Ntotalmol];
	int *movestat = new int[Ntotalmol];
	double h;

	/*The number of species does not include all separate interfaces, because it is for diffusible species*/

	/*Iterate over time steps until you hit a max time*/
	int mu;
	int flag2;
	double prob;
	double sum;
	double xchg, ychg, zchg;
	int icom;
	int pro_type, mp;
	int whichspecie;
	double rerand;
	double hfact;
	double dx, dy, dz;
	double dist2;
	double box_x = plist.xboxl;
	double box_y = plist.yboxl;
	double box_z = plist.zboxl; //nm
	double xtot, ytot, ztot;
	int nfree;
	int p, i1, i2;
	int np;
	double r2, r;
	double maxsep2 = bindrad[0] * bindrad[0]; //plist.maxsep2;
	cout << "squared distance cutoff: " << maxsep2 << endl;
	int iind, iind2, ppart;
	
	int statwrite = plist.statwrite;
	int configwrite = plist.configwrite;
	int restartwrite = plist.grwrite;
	// we have to do the following two lines so that timestat file and restart file are in time-sync
	if(restartwrite % statwrite != 0){
	  int tmp=int(restartwrite/statwrite);
	  if(tmp==0)tmp=1;
	  restartwrite=tmp*statwrite;
	}

	int *nBoundPairs=new int[Nprotypes*Nprotypes];//Array to hold how many protein pairs have bound together.
	for(i=0;i<Nprotypes*Nprotypes;i++)
	  nBoundPairs[i]=0;

	sprintf(fnameComplexXYZ, "ComplexComponents_Np%d_Ni%d.dat", plist.Nprotypes, Nifaces);
	ofstream compout(fnameComplexXYZ);
	sprintf(fnameProXYZ, "ProteinW_COM_Np%d_Ni%d.xyz", plist.Nprotypes, Nifaces);
	ofstream proout(fnameProXYZ);
	sprintf(fnameProXYZ, "debug_dump.dat");
	ofstream tfile(fnameProXYZ);
	sprintf(fnameProXYZ, "ComplexHistogram_Np%d.Ni%d.dat", plist.Nprotypes, Nifaces);
	ofstream assemblyfile(fnameProXYZ);
	sprintf(fnameProXYZ, "MonomerDimerFile.dat");
	ofstream dimerfile(fnameProXYZ);
	sprintf(fnameProXYZ, "MeanComplexSize.dat");
	ofstream meanfile(fnameProXYZ);

	/*Establish the largest Rmax between binding pairs, so the cells are big enough
	 Loop over all reactions and include the distance to travel, also include the distance from interfaces
	 to the protein COM, because the bin is only defined by the protein COM!!!*/
	
	double Rmaxlimit=0;
	double Rmax_radius=0;//how much of the Rmax is due to distance from COM to binding interfaces

	Rmaxlimit=set_Rmaxlimit( plist,  wholep,  bindrad,  i_home,  p_home,  Rlist, rxtype,  coordcont, Rmax_radius);
	
	double Rmax_add=Rmaxlimit-Rmax_radius;
	double deltat = plist.deltat;
	double cellSizeLength=Rmaxlimit*2.0;//Instead of using the shortest length reachable by diffusion, make it a bit larger. One reason is for large complexes, Deff can becomre larger and cause larger translations than anticipated at the beginning. 
	double tx, ty, tz;

	/*cells*/
	double *M = new double[9];
	int myc1, myc2;
	int maxnbor = 13; //geometry of cube
	int Nx, Ny, Nz;
	
	Nx = int(floor(plist.xboxl / cellSizeLength));//-1;//minus 1 is somewhat arbitrary, should be safe without it.
	Ny = int(floor(plist.yboxl / cellSizeLength));//-1;
	int scale_down=1;//have not yet optimized speed for number of cells.
	if(Rmaxlimit/Rmax_add>2)scale_down=1;
	else scale_down=4;
	Nz = max(2,int(floor(plist.zboxl / cellSizeLength))/scale_down );
	//	Nz=400;
	double cellx = box_x / (Nx * 1.0);
	double celly = box_y / (Ny * 1.0);
	double cellz = box_z / (Nz * 1.0);
	int Ncell = Nx * Ny * Nz;
	if (cellx < cellSizeLength){ 
	  cout << "CELL SIZE IS TOO SMALL " << " Rmax: "<<cellSizeLength<<" cellsize: "<<cellx<<endl;
	  Nx = max(2, int(floor(plist.xboxl / cellSizeLength))/2);
	}
	if (celly < cellSizeLength){ 
	  Ny = max(2, int(floor(plist.yboxl / cellSizeLength))/2); 
	}
	if (cellz < cellSizeLength){ 
	  Nz = max(2,int(floor(plist.zboxl / cellSizeLength))/2 ); 
	}
	Ncell = Nx * Ny * Nz;
	cellx = box_x / (Nx * 1.0);
	celly = box_y / (Ny * 1.0);
	cellz = box_z / (Nz * 1.0);	
	cout << "Nx: " << Nx << " Ny " << Ny << " Nz " << Nz << " Ncell " << Ncell << endl;
	double div=2.0;
	double maxPairs=0.5*Ntotalmol*Ntotalmol;
	while(Ncell > maxPairs){//Ncell*maxnbor
	 	cout<< "Maximum number of cells " << Ntotalmol*Ntotalmol*0.5 << " exceeded!"<< endl;
	 	cout<< "Scaling down number of cells" <<endl;
		
		Nx = max(2, int(floor(plist.xboxl / cellSizeLength)/div));
		Ny = max(2, int(floor(plist.yboxl / cellSizeLength)/div));
		Nz = max(2,int(floor(plist.zboxl / cellSizeLength)/(2.0*div)) );
	
	 	Ncell = Nx * Ny * Nz;
	 	
	 	cellx = box_x / (Nx * 1.0);
	 	celly = box_y / (Ny * 1.0);
	 	cellz = box_z / (Nz * 1.0);
		div=div+1;
	 }
	
	cout << "Nx: " << Nx << " Ny " << Ny << " Nz " << Nz << " New Ncell " << Ncell << endl;
	cout << "Rmaxlimit: " << Rmaxlimit <<" Set Min CellSizeLength: "<<cellSizeLength<<  " final cell size length, x: " << cellx << endl;
	
	
	/*Partition full cell volume into cells, to speed up search for binding partners*/
	/*Each cell has 26 neighbors (unless at boundary and reflecting BC is used).
	  Since we will loop over all cell pairs, only keep track of half your neighbors.
	 */
	
	int *Nnbor=new int[Ncell];
	int *nbor = new int[Ncell * maxnbor];
	//	int *nborrev = new int[Ncell * maxnbor];
	int *npb = new int[Ncell];

	vector< vector<int> > binlist(Ncell);
	//int *binlist=new int[Ncell*MAXPERBIN];
	
	vector<int> disslist(10); //holds index of dissociated proteins. Initialize size to some value.
	int ndiss;
	int mybin;
	int mybinind;
	int c;
	cell_neighbor_list(Nx, Ny, Nz, maxnbor, nbor, Nnbor);
	cout << "N cell pairs (max is with PBC): " << Ncell * maxnbor << endl;

	
	

	/***************************/

	cout << "deltat: " << plist.deltat << endl;
	long int it, iterass = 0;
	double pnormval, pirrval;
	ifstream restartf;
	plist.ntotalcomplex = Ntotalmol;
	if (plist.restart == 1) {
		/*update status of each protein and complex.*/

		restartf.open(argv[3]);
		iterass = read_restartCELL(restartf, Ntotalmol, bases, plist, ind_com, Ncopy, wholep);
		restartf.close();
		/*get complex com, complex radius, and complex diffusion*/
		update_complex_all(plist.ntotalcomplex, ind_com, bases);
	}
	long int Nit = int(plist.Nit);
	if (Nit > 2.147E9) {
		cout << "ITERATIONS EXCEEDS INTEGER MAXIMUM! Exiting..." << endl;
		exit(1);
	}
	int s1;
	cout << "Ntotal complexes: " << plist.ntotalcomplex << endl;
//	write_complex(compout, plist, ind_com, 0);
//	write_protein_iface_short(proout, plist, bases, Ncopy, 0, wholep, names);
//	write_dcd(proout, plist, bases, Ncopy, 0, wholep, infofilenames);

	int amol, df;
	double us_to_s = 1E-6;


	double **traj = new double *[Ntotalmol];
	double **trajR = new double *[Ntotalmol];
	for (i = 0; i < Ntotalmol; i++) {
		traj[i] = new double[3];
		trajR[i] = new double[3];
	}

	double R1;

	double r0, passoc;
	int veclen;

	int maxnbort;


	/*Below use einstein-stokes to define D based on the particle
	 radius and the Temperature and viscosity.
	 To enforce the D read in from file, calculated via, e.g. bead models,
	 set scale>1 below.
	 */
	double Temp = 293; //K
	double nu = 0.001; //kg/(m*s)
	double scale = 3.0; //greater the one to correct for non-spherical
	double crad = wholep[0].radR;
	if (crad == 0)
		crad = 1;
	cout << "clathrin radius: " << crad << " nm. " << endl;
	plist.pretrans = trans_prefactor(Temp, nu, scale, crad, wholep[0].Dx);
	plist.prerot = rot_prefactor(Temp, nu, scale, crad, wholep[0].Drx);
	cout << "Diffusion prefactors: " << plist.pretrans << ' ' << plist.prerot << endl;
	
	
	int loop = 0;
	
	int *Nsum = new int[Nprotypes];
	Nsum[0] = 0;
	for (i = 1; i < Nprotypes; i++) {
		Nsum[i] = Nsum[i - 1] + Ncopy[i - 1];
	}
	//	int size = Ncopy[0] * Ncopy[1];
	cout << "Ncopy[0]: " << Ncopy[0] << endl;

	int *nprevpart = new int[Ntotalmol];
	int *ncurrpart = new int[Ntotalmol];
	int **prevlist = new int*[Ntotalmol];
	int **currlist = new int*[Ntotalmol];
	int **prevmyface = new int*[Ntotalmol];
	int **currmyface = new int*[Ntotalmol];
	int **prevpface = new int*[Ntotalmol];
	int **currpface = new int*[Ntotalmol];
	double **prevnorm = new double*[Ntotalmol];
	//int **previter=new int*[Ntotalmol];
	double **ps_prev = new double*[Ntotalmol];
	double **prevsep = new double*[Ntotalmol];
	double **currprevnorm = new double*[Ntotalmol];
	//int **previter=new int*[Ntotalmol];
	double **currps_prev = new double*[Ntotalmol];
	double **currprevsep = new double*[Ntotalmol];

	/*This size MAXNORM is keep track of how many particles are in each proteins reaction
	 zone at one time. Should be set large to buffer for fluctuations to large number
	 of partners, even though each interface (and protein) should ideally have only 1
	 partner in its reaction zone at each step.
	 */
	int MAXNORM = 200;
	int s;
	int ssave;
	for (i = 0; i < Ntotalmol; i++) {
		nprevpart[i] = 0;
		ncurrpart[i] = 0;

		currlist[i] = new int[MAXNORM];
		prevlist[i] = new int[MAXNORM];
		currmyface[i] = new int[MAXNORM];
		prevmyface[i] = new int[MAXNORM];
		currpface[i] = new int[MAXNORM];
		prevpface[i] = new int[MAXNORM];
		prevnorm[i] = new double[MAXNORM];
		ps_prev[i] = new double[MAXNORM];
		prevsep[i] = new double[MAXNORM];
		currprevnorm[i] = new double[MAXNORM];
		currps_prev[i] = new double[MAXNORM];
		currprevsep[i] = new double[MAXNORM];
		//previter[i]=new int[MAXNORM];
	}
	//cerr<<"allocated mem "<<endl;
	int myface, pface;
	int place2;
	double rtol = 1E-10;
	int flag;
	double tol=1E-20;
	int p1, p2;
	double probvec1, p0_ratio, currnorm;

	cout << "Set prevnorms to one " << endl;
	for (i = 0; i < Ntotalmol; i++) { //previously this was Ncopy[0] by Dr. Johnson
		currlist[i][0] = 0;
		for (j = 0; j < MAXNORM; j++) {

			prevnorm[i][j] = 1.0;

			//previter[i][j]=-1;
			ps_prev[i][j] = 0;
			prevsep[i][j] = 0;
		}
	}
	
 
	double Nacurr;
	double pact;
	double tmpx, tmpy, tmpz;

	int pnew;
	int c1, c2;
	int ci1, ci2;
	double tremain, tevent;
	int mu_ret;
	int rxn1;
	int go;
	double rate;
	int cancel;
	double meanComplexSize;

	int nc1, nc2;
	double sep, ratio;
	double aexp;
	double bexp;

	int flagsep;
	
	int pp, qq, nb, hh;

	double maxrad = 150;
	int radbins = 3000;
	double delrad = maxrad / (1.0 * radbins);
	int ind_rad;
	double rad2, rad;

	int Ncsave = plist.ntotalcomplex;
	int Nccurr;

	double x0;
	int proa, pro2;
	double currx, curry, currz;
	int Nrep = 1;
	int rep = 0;
	int totpercell = 0;
	int dub;
	
	/*Determine which bin the particles are in*/
	// for(i=0;i<Ntotalmol;i++){
	//   cout <<"Protein : "<<i<<" Dz: "<<bases[i].Dz<<" radius: " <<bases[i].radR<<" complex: "<<ind_com[i].Dz<<" radius: "<<ind_com[i].radR<<endl;
	// }
	get_bin2(plist, bases, cellx, celly, cellz, Nx, Ny, Nz, binlist, npb, Ncell, ind_com, 0);
	for(i=0;i<Ncell;i++){
	  //cout << "cell: " << i << " npercell: " << npb[i] << endl;
	  
	  totpercell += npb[i];
	}
	cout << "Done defining bins. Total mols across cells " << totpercell << endl;
	std::vector<int> proPairlist;
	sprintf(fnameProXYZ, "NboundPairFile.dat");
	ofstream pairOutfile(fnameProXYZ);
	init_NboundPairs(nBoundPairs,  proPairlist, pairOutfile,  plist, infofilenames, wholep);
	it=0;
	write_NboundPairs(nBoundPairs, proPairlist,  pairOutfile, 0, plist);
	
	meanComplexSize=calc_complex_hist(plist.ntotalcomplex, ind_com, bases, assemblyfile, 0, plist, infofilenames, Ncopy);
	meanfile <<it*plist.deltat<<'\t'<<meanComplexSize<<endl;
	init_print_dimers(dimerfile, 0, plist, infofilenames);//Initialize header to dimerfile at it=0
	print_dimers(plist.ntotalcomplex, ind_com, bases, dimerfile, it, plist, infofilenames, Ncopy);
	cerr<<" Done writing out to new output files "<<endl;
	
	//calc_complex_hist(plist.ntotalcomplex, ind_com, bases, assemblyfile, 0, plist, infofilenames, Ncopy);
	/*Write out initial coordinate files.*/
	if (iterass == 0) {
	  ofstream restart(tname);
	  write_restart(restart, wholep, bases, plist, Ncopy, it + iterass, ind_com, deltat);
      //	  write_timestat(timestatfile, wholep, bases, plist, Ncopy, iterass, ind_com, deltat, Nprotypes); ///was assemblyfile
	  //write_timestat2(timestatfile, molectypesfile, timestatfiletext, wholep, bases, plist, Ncopy, it + iterass, ind_com, deltat, Nprotypes, molectypes, currentnumberofmolectypes);
	  //Have to uncomment these when you want to deal with restart files (coords data)
	  
	  //		write_protein_iface_short(proout, plist, bases, Ncopy, (it + iterass), wholep, names);
	  //  write_dcd(proout, plist, bases, Ncopy, 0, wholep, infofilenames);
	  cerr <<"write complex components: "<<endl;
	  write_complex_components(plist.ntotalcomplex, ind_com, bases, compout, 0, plist, infofilenames);
	  write_protein_iface(proout, plist, bases, Ncopy, 0, wholep);
	  restart.close();
	}
	

	
	
	
	
	// cout <<"Initial coordinates: "<<endl;
	// for(i=0;i<Ntotalmol;i++){
	//   write_crds(bases, i);
	//   cout <<"Protein: "<<i<<" Ninterface: "<<bases[i].ninterface<<" x cord: "<<bases[i].x[0]<<endl;
	  
	// }
	/*********************************/
	/*BEGIN SIMULATION*/
	/******************************/
	
	for (it = 1; it < Nit + 1; it++) {
	  
	  for (i = 0; i < Ntotalmol; i++) {
	    ncross[i] = 0;
	    ncrosscom[i]=0;
	    movestat[i] = 0;
	  }
	  disslist.clear();
	  /***************/
	  /*Test dissociation*/
	  
	  dissociate_sigma_com( Ntotalmol,  bases,  ind_com,  plist,  p_home,  i_home,  Rlist,  Nmyrxn, myrxn,  bindrad,  Ncoup,  mycoupled,  irandmax,  disslist, ncross, kr, it, ncrosscom, movestat, nBoundPairs);
	 
   	  /*Produce Checkpoint files*/
	  if (it % (restartwrite) == 0) { 
		  // precautionary: need to have this otherwise the restart file gets garbled when we're saving for a very large number of molecules and simulation abruptly were terminated
	    std::ifstream src0("restartsW.out", std::ios::binary);
	    std::ofstream dst0("SAVEDrestartsW.out", std::ios::binary);
	    dst0 << src0.rdbuf();
	    std::ifstream src1("timestat.dat", std::ios::binary);
	    std::ofstream dst1("SAVEDtimestat.dat", std::ios::binary);
	    dst1 << src1.rdbuf();

//	    Have to uncomment these when you want to deal with restart files (coords data)
//	    std::ifstream src2(fnameComplexXYZ, std::ios::binary);
//	    std::ofstream dst2("SAVEDfnameComplexXYZ.dat", std::ios::binary);
//	    dst2 << src2.rdbuf();
//	    std::ifstream src3(fnameProXYZ, std::ios::binary);
//	    std::ofstream dst3("SAVEDfnameComplexXYZ.dat", std::ios::binary);
//	    dst3 << src3.rdbuf();
	  }  

	  
	  get_bin2(plist, bases, cellx, celly, cellz, Nx, Ny, Nz, binlist, npb, Ncell, ind_com, it);
	  // totpercell=0;
	  // for(i=0;i<Ncell;i++)
	  //   totpercell += npb[i];
	  // cout <<"Tot per bin check: "<<totpercell<<" Bin 208: "<<bases[208].mybin<<" bin 297: "<<bases[297].mybin<<endl;
	  /*Either remove dissociated proteins from the binlist
	    so they do not try to react again this time step, or
	    ensure that they do not move again (D=0) and they
	    will not react, the other proteins will just avoid
	    overlapping with them.
	  */
	  for (i = 0; i < disslist.size(); i++) {
	    p1 = disslist[i];
	    k = bases[p1].mycomplex;
	    cout <<" In disslist: p1: "<<p1<<endl;
	    traj[k][0]=0;
	    traj[k][1]=0;
	    traj[k][2]=0;
	    /*Comment out below to allow dissociated proteins to avoid overlap, but
	     need to update evaluate binding, to test if movestat!=2*/
	    /*
	    mybin = bases[p1].mybin;
	    mybinind = bases[p1].mybinind;
	    
	    pnew = binlist[mybin * MAXPERBIN + npb[mybin] - 1];
	    bases[pnew].mybinind = mybinind;
	    binlist[mybin * MAXPERBIN + mybinind] = pnew;
	    npb[mybin] -= 1;
	    */
	  }
	  
	  /*Keep track of partners in reaction zone for reweighting*/
	  for (i = 0; i < Ntotalmol; i++) {
	    ncurrpart[i] = 0;
	  }
	  
	  /*Measure separations between proteins in neighboring cells to identify
	    all possible reactions.
	  */
	  
	  for (c = 0; c < Ncell; c++) {
	    //			cout<<npb[c]<<endl;
	    for (pp = 0; pp < npb[c]; pp++) {
	      i = binlist[c][pp];
	      
	      /*Test bimolecular reactions!*/
	      nfree = bases[i].nfree;
	      	      
	      if (nfree > 0) {
		/*first loop over proteins in your same cell.*/
		for (qq = pp + 1; qq < npb[c]; qq++) {
		  j = binlist[c][qq];
		  
		  
		  evaluate_binding_pair_com(i, j, bases,  ind_com,  plist,  wholep,  DDtableindex,  numpartners,  Speclist,  i_home,  ncross,  probvec,  TBLID, bindrad,  nprevpart,  prevlist,  prevmyface, prevpface,  prevnorm,  ncurrpart,  currlist,  currmyface,  currpface,  ps_prev, prevsep,  currprevnorm,  currps_prev,  currprevsep, myrxn, traj, crosspart, crossint, cross_rxn, kr, it, contsur, contnorm, contpir, MAXALWTBL, ncrosscom, movestat);
		  
		} //loop over protein partners in your same cell
		
		/*Now loop over all neighboring cells, and all proteins in those cells.
		  for PBC, all cells have maxnbor neighbor cells. For reflecting, edge have fewer.
		*/
		
		
		for (hh = 0; hh < Nnbor[c]; hh++) {
		  nb = nbor[c * maxnbor + hh];
		  for (qq = 0; qq < npb[nb]; qq++) {
		    j = binlist[nb][qq];

		    evaluate_binding_pair_com(i, j, bases,  ind_com,  plist,  wholep,  DDtableindex,  numpartners,  Speclist,  i_home,  ncross,  probvec,  TBLID, bindrad,  nprevpart,  prevlist,  prevmyface, prevpface,  prevnorm,  ncurrpart,  currlist,  currmyface,  currpface,  ps_prev, prevsep,  currprevnorm,  currps_prev,  currprevsep, myrxn, traj, crosspart, crossint, cross_rxn, kr, it, contsur, contnorm, contpir, MAXALWTBL, ncrosscom, movestat);
		    
		    
		  } //loop over all proteins in this neighbor cell
		} //loop over all neighbor cells
	      } //if protein i is free to bind
	    } //loop over all proteins in initial cell
	  } //End looping over all cells.
	  
	  /*Now that separations and reaction probabilities are calculated,
	    decide whether to perform reactions for each protein.
	  */
	  
	  for (i = 0; i < Ntotalmol; i++) {
	    
	    /*Skip any proteins that just dissociated during this time step*/

	    if (ncross[i] > 0) {
	      
	      	      	      
	      /*Evaluate whether to perform a reaction with protein i, and with whom. Flag=1 means
		reaction is performed. Returns correct ci1 and ci2 for this rxn.
	      */
	      
	      /*Loop over all reactions individually, instead of summing probabilities*/
	      flag = choose_one_reaction_loop( i, ncross, crosspart, probvec, ci1, ci2, cross_rxn, crossint, irandmax);


	      /*FOR DEBUGGING, FORCE ASSOCIATION TO OCCUR!!*/

	      /*flag = 0; //we chose an association reaction
	      ci1 = 0;
	      p1=i;
	      p2 = crosspart[p1][ci1];
	      double pmatch = probvec[p1][ci1];
	      if(pmatch>0)flag=1; 
	      for (j = 0; j < ncross[p2]; j++) {
		if (probvec[p2][j] == pmatch) {
		
		  ci2 = j;
		  if (cross_rxn[p1][ci1] == cross_rxn[p2][ci2])
		    j = ncross[p2]; //break from loop, otherwise the next prob 0 point could be your reactant
		}
	      }
	      */
	    	    
	      if (flag == 1) {
		
		perform_association_sigma_com(i,  bases,  ind_com, plist,  crosspart,  cross_rxn, crossint,  Rlist,  i_home,  probvec,  Ncoup,  mycoupled,  p_home,  ci1,  ci2,  it, movestat, ncross, bindrad[cross_rxn[i][ci1]], ncrosscom, traj,trajR, nBoundPairs);
		
		
	      } else {
		
		/*This protein will not associate in this time step.
		  For this Sweeping version, just move the particle, don't check for overlap until everyone is done.
		  Store new move in the traj vector so you still know the original position store in bases[]
		*/
		
		/*movestat of zero means no traj value is selected.
		  movestat=1 means traj is selected, but particles have not moved
		  movestat=2 means particles have moved
		*/
		if (movestat[i] == 0) {
		  k = bases[i].mycomplex;
		  dx = sqrt(2.0 * deltat * ind_com[k].Dx) * GaussV();
		  dy = sqrt(2.0 * deltat * ind_com[k].Dy) * GaussV();
		  dz = sqrt(2.0 * deltat * ind_com[k].Dz) * GaussV();
		  
		  traj[k][0] = dx;
		  traj[k][1] = dy;
		  traj[k][2] = dz;
		  
		  tx = sqrt(2.0 * deltat * ind_com[k].Drx) * GaussV();
		  ty = sqrt(2.0 * deltat * ind_com[k].Dry) * GaussV();
		  tz = sqrt(2.0 * deltat * ind_com[k].Drz) * GaussV();
		  
		  trajR[k][0] = tx;
		  trajR[k][1] = ty;
		  trajR[k][2] = tz;
		  
		  rotationEuler(trajR[k][0], trajR[k][1], trajR[k][2], M);
		  reflect_traj_complex_rad_rot(i, bases, ind_com, plist, traj[k], trajR[k], M);
		  /*this displacement will apply to all the proteinss in this complex k.*/
		  for (j = 0; j < ind_com[k].mysize; j++) {
		    mp = ind_com[k].plist[j];
		    movestat[mp] = 1;
		  }
		  
		}
		
		/*Set probability of this protein to zero in all reactions so it doesn't try to
		  react again but the partners still will avoid overlapping.
		*/
		remove_one_prob_all(i, ncross, crosspart, probvec, cross_rxn, crossint);
		
	      }
	      
	    } else if (ncross[i] > -1) {
	      
	      /*this protein has ncross=0,
		meaning it neither dissociated nor tried to associate.
		however, it could have movestat=2 if it is part of a multi-protein
		complex that already displaced.
	      */
	      if (movestat[i] == 0) {
		/*don't move this protein if it already moved (movestat=2), or attached to someone who will
		  potentially have to resample their position (movestat=1)*/
		k = bases[i].mycomplex;
		dx = sqrt(2.0 * deltat * ind_com[k].Dx) * GaussV();
		dy = sqrt(2.0 * deltat * ind_com[k].Dy) * GaussV();
		dz = sqrt(2.0 * deltat * ind_com[k].Dz) * GaussV();
		
		traj[k][0] = dx;
		traj[k][1] = dy;
		traj[k][2] = dz;
		
		tx = sqrt(2.0 * deltat * ind_com[k].Drx) * GaussV();
		ty = sqrt(2.0 * deltat * ind_com[k].Dry) * GaussV();
		tz = sqrt(2.0 * deltat * ind_com[k].Drz) * GaussV();
		
		trajR[k][0] = tx;
		trajR[k][1] = ty;
		trajR[k][2] = tz;
		
		/*Don't automatically move the protein, in case another part of your protein complex has overlap
		 */
		rotationEuler(trajR[k][0], trajR[k][1], trajR[k][2], M);
		reflect_traj_complex_rad_rot(i, bases, ind_com, plist, traj[k], trajR[k], M);
		  
		
		//if other proteins in this complex have already moved, don't need to sample them as well;
		for(j=0;j<ind_com[k].mysize;j++){
		  mp=ind_com[k].plist[j];
		  movestat[mp]=1;
		}
	    
	      }
	      
	    }
	    
	  } //done testing all molecules for bimolecular reactions
	  /*Now we have to check for overlap!!!*/
	  for (i = 0; i < Ntotalmol; i++) {
	    
	    /*
	      Now track each complex (ncrosscom), and test for overlap of all proteins in that complex before
	      performing final position updates. 
	    */
	    
	    if (ncrosscom[bases[i].mycomplex] > 0) {
	      if(movestat[i]!=2){
		/*For any protein that overlapped and did not react, check whether it overlaps with its partners,
		  do all proteins in the same complex at the same time.
		  Also, if both proteins are stuck to membrane, only do xy displacement, ignore z
		*/
		if(ind_com[bases[i].mycomplex].Dz<tol){
		  //cout <<"Protein memtest: "<<i<<" ncross: "<<ncross[i]<<" ncrosscom: "<<ncrosscom[bases[i].mycomplex]<<endl;
		  sweep_separation_complex_rot_memtest(deltat, i, bases, ind_com, ncross, crosspart, crossint, cross_rxn, traj, probvec, plist, movestat, i_home, bindrad, trajR, M, Rlist, it);
		}else
		  sweep_separation_complex_rot(deltat, i, bases, ind_com, ncross, crosspart, crossint, cross_rxn, traj, probvec, plist, movestat, i_home, bindrad, trajR, M, Rlist, it);
		
		
	      }
	    } else if (ncrosscom[bases[i].mycomplex] == 0) {
	      
	      if (movestat[i] != 2) {
		
		/*For proteins with ncross=0, they either moved independently, or their displacements
		  were selected based on the complex they were part of, and they may not yet been moved.
		*/
		/*Update their trajectories first to fit in the box. Then move the particles.*/
		k=bases[i].mycomplex;
		rotationEuler(trajR[k][0], trajR[k][1], trajR[k][2], M);
		reflect_traj_complex_rad_rot(i, bases, ind_com, plist, traj[k],trajR[k], M);
		move_rot_proteins(i, bases, ind_com, traj[k], movestat, trajR[k], M);
		//reflect_complex_rad_rot(i, bases, ind_com, plist.xboxl, plist.yboxl, plist.zboxl);
		
	      }
	      
	    }
	  }
	  
	  if (it % statwrite == 0) {
	    cout <<endl;
	    cout <<" timestep: "<<it*deltat<<" iteration: "<<it<<" total complexes: "<<plist.ntotalcomplex<<endl;
	    cout <<endl;
	    meanComplexSize=calc_complex_hist(plist.ntotalcomplex, ind_com, bases, assemblyfile, it, plist, infofilenames, Ncopy);
	    print_dimers(plist.ntotalcomplex, ind_com, bases, dimerfile, it, plist, infofilenames, Ncopy);
	    write_NboundPairs(nBoundPairs, proPairlist,  pairOutfile, it, plist);
	    meanfile <<it<<'\t'<<meanComplexSize<<endl;
	    	    
	    if(it % configwrite ==0)
	      write_protein_iface(proout, plist, bases, Ncopy, it, wholep);
	  }
	  
	  if (it % restartwrite == 0) {//restart file written here
	    ofstream restart(tname);
	    write_restart(restart, wholep, bases, plist, Ncopy, it + iterass, ind_com, deltat);
	    restart.close();
	    //Have to uncomment these when you want to deal with restart files (coords data)
	    //		write_complex(compout, plist, ind_com, deltat*(it + iterass));
	    //		write_protein_iface_short(proout, plist, bases, Ncopy, (it + iterass), wholep, names);
	  }
	  
	  /*Now replace all currsep as prevsep, to keep track of
	    reweighting values for next time step*/
	  
	  for (i = 0; i < Ntotalmol; i++) {
	    nprevpart[i] = ncurrpart[i];
	    for (s = 0; s < ncurrpart[i]; s++) {
	      prevlist[i][s] = currlist[i][s];
	      prevmyface[i][s] = currmyface[i][s];
	      prevpface[i][s] = currpface[i][s];
	      prevnorm[i][s] = currprevnorm[i][s];
	      ps_prev[i][s] = currps_prev[i][s];
	      prevsep[i][s] = currprevsep[i][s];	      
	    }
	  }	  
  
	} //end iterating over time steps
	
	stop_timer(&totaltime);
	cout << timer_duration(totaltime) << " total time " << endl;
	
	/*Write out final result*/
	cout << "End Main, complete run " << endl;

} //end main
