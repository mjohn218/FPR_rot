#ifndef __INTEGRATE_RD_H
#define __INTEGRATE_RD_H

#include "rd_class_defs.h"
//#include "cluster.h"

/*define ka for membrane binding that recovers kon_3D*/
double ka_rec(double s, double ka3D, double Dasol, double Dbsol, double fourpi) ;

/*Determine which cell in the partitioned volume the particle is in*/
void get_bin2(Parms plist, Fullmol *bases, double cellx, double celly, double cellz, int Nx, int Ny, int Nz, std::vector<std::vector <int> > &binlist, int *npb, int Ncell, Complex *ind_com, int iter);
void put_back_inside(int p1, Parms plist, Fullmol *bases, int *npb, int Ncell, Complex *ind_com, int iter, int &loopit);

void ayudame(int p1, int p2, Fullmol *bases);

/*Choose new position for particles, avoid all other partners in reaction zone*/
void sweep_separation(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad);
void sweep_separation_all(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, int **Rlist);
void sweep_separation_allPBC(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, int **Rlist);
void sweep_separation_rot_allPBC(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, double **trajR, double *M, int **Rlist);
void sweep_separation_rot_all(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, double **trajR, double *M, int **Rlist);
void sweep_separation_complexPBC(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, int **Rlist);
void sweep_separation_complex_rot_PBCCELL(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, double **trajR, double *M, int **Rlist, int iter);
void sweep_separation_complex_rot(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, double **trajR, double *M, int **Rlist, int iter);
void sweep_separation_complex_rotPBC(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, double **trajR, double *M, int **Rlist, int iter);
void sweep_separation_complex_rot_memtest_PBCCELL(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, double **trajR, double *M, int **Rlist);
void sweep_separation_complex_rot_memtest(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, double **trajR, double *M, int **Rlist, int iter);
void sweep_separation_complex_rot_memtestPBC(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, double **trajR, double *M, int **Rlist, int iter);
void cluster_sweep_separation_complex_rot_memtestPBC(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, double **trajR, double *M, int **Rlist, double *kr, int iter) ;
void sweep_separation_rot_allPBCCELL(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, double **trajR, double *M, int **Rlist);
void sweep_separation_rot_com_allPBC(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, double **trajR, double *M, int **Rlist);
void sweep_separation_vr_allPBC(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, int **Rlist, double **Ftraj);
void islander(int *movestat, int *islandsize, int *islandvec, int newsizeoverlap, int Ntotalmol, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **neighbors, int &maxislandsize, double &meanislandsize, int &islandnum);
void islander_stats(int *islandsize, int *islandvec, int islandnum, int Ntotalmol, int &maxislandsize, double &meanislandsize);
int isneighbor(int p1, int p2, int **neighbors);
int check_overlap(double &dr2, int p1, int p2, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, Parms plist, int *ihome, double *bindrad, double **trajR, double *M, int **Rlist);
int check_overlapPBC(double &dr2, int p1, int p2, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, Parms plist, int *ihome, double *bindrad, double **trajR, double *M, int **Rlist);
void all_sweep_separation_complex_rot_memtest_PBCCELL(int &overlaps, int myislandsize, int **neighbors, int islandid, int Ntotalmol, int *islandvec, double deltat, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, double **trajR, double *M, int **Rlist);
void all_sweep_separation_complex_rot_memtest_PBC(int &overlaps, int myislandsize, int **neighbors, int islandid, int Ntotalmol, int *islandvec, double deltat, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, double **trajR, double *M, int **Rlist);

void cluster_one_complex(int k1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, Parms plist, double *bindrad,  int **Rlist, std::vector<ClusterPair> &pairList, std::vector<int> &finished, double *kr) ;
void define_cluster_pairs(int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn,  Parms plist,  double *bindrad,  int **Rlist, std::vector<ClusterPair> &pairList, double *kr) ;
void resample_traj(int currStop, std::vector<ClusterPair> &pairList, double **traj, double **trajR, Complex *ind_com, double deltat, double *M, int *movestat);
void cluster_priority_onlyPBC(std::vector<ClusterPair> &origPairList, double deltat,  Fullmol *bases, Complex *ind_com, int *ncross,  double **traj,  Parms plist, int *movestat, int *ihome, double **trajR, double *M,  int iter) ;

/*get distances from GF*/
void sampleGF_separation_allPBC(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, int **Rlist, double **rsample);
void stretch_multiple_cross(int p1, int *ncross, Complex *ind_com, Fullmol *bases, double **traj, Parms plist, int **cross_part, double **rsample, int *movestat);
void stretch_pro_com_coordsPBC(int p1, int p2, Complex *ind_com, Fullmol *bases,Parms &plist, double Rnew);
void stretch_p2_com_coordsPBC(int p1, int p2, Complex *ind_com, Fullmol *bases,Parms &plist, double Rnew, int **cross_part);


/*Don't move again*/
void set_movestat_zero(int p1, Fullmol *bases, Complex *ind_com, int *movestat);

/*Recalculate COMs for the complexes, and for some the D and radius*/
void update_complex_all(int Ncomplex, Complex *ind_com, Fullmol *bases);
void update_one_com_only(int c1, Complex *ind_com, Fullmol *bases);
void update_one_com_onlyPBC(int c1, Complex *ind_com, Fullmol *bases, Parms plist,int flagbndry);
void update_one_com_onlyPBCCELL(int c1, Complex *ind_com, Fullmol *bases, Parms plist, int flagbndry);

/*If a protein has reacted, delete it from the list of other possible partners*/
void remove_protein_reaction(int p1, int *ncross, int **crosspart, double **prob,int **cross_rxn, int **cross_int, int *ncrosscom, Fullmol *bases);
//void remove_reaction_all(int p1, int p2, int *ncross, int **crosspart, double **prob, int **cross_rxn, int **cross_int, int *ncrosscom, Fullmol *bases);
void remove_reaction_all(int p1, int p2, int *ncross, int **crosspart, double **prob, int ci1, int ci2, int **cross_rxn, int **cross_int);
void remove_reaction_all_ncom(int p1, int p2, int *ncross, int **crosspart, double **prob, int **cross_rxn, int **cross_int, int *ncrosscom, Fullmol *bases);
void remove_one_prob_all(int p1, int *ncross, int **crosspart, double **prob, int **cross_rxn, int **cross_int);

void remove_reaction_all_GFr(int p1, int p2, int *ncross, int **crosspart, double **prob, int **cross_rxn, int **cross_int, double **cross_sep, double **rsample, int *ncrosscom, Fullmol *bases);

/*Measure separation between a particle pair*/
int get_distance(Fullmol *bases, double **traj, int p1, int p2, double deltat, double bindrad, int *ncross, int **crosspart, int **crossint, int i1, int i2, int mu, int **cross_rxn, int *ihome, double iter, double Rmax, double &sep, double &R1, int *ncrosscom);
int get_distancePBC(Fullmol *bases, double **traj, int p1, int p2, double deltat, double bindrad, int *ncross, int **crosspart, int **crossint, int i1, int i2, int mu, int **cross_rxn, int *ihome, double iter, double Rmax, double &sep, double &R1, int *ncrosscom, Parms &plist);
int get_distancePBC_2D(Fullmol *bases, double **traj, int p1, int p2, double deltat, double bindrad, int *ncross, int **crosspart, int **crossint, int i1, int i2, int mu, int **cross_rxn, int *ihome, double iter, double Rmax, double &sep, double &R1, int *ncrosscom, Parms &plist);
int get_distance_comPBC(Fullmol *bases, double **traj, int p1, int p2, double deltat, double bindrad, int *ncross, int **crosspart, int **crossint, int i1, int i2, int mu, int **cross_rxn, int *ihome, double iter, double Rmax, double &sep, double &R1, Parms &plist);
int get_distancePBCCELL(Fullmol *bases, double **traj, int p1, int p2, double deltat, double bindrad, int *ncross, int **crosspart, int **crossint, int i1, int i2, int mu, int **cross_rxn, int *ihome, double iter, double Rmax, double &sep, double &R1, Parms &plist, int *ncrosscom);
int get_distancePBCCELL_2D(Fullmol *bases, double **traj, int p1, int p2, double deltat, double bindrad, int *ncross, int **crosspart, int **crossint, int i1, int i2, int mu, int **cross_rxn, int *ihome, double iter, double Rmax, double &sep, double &R1, Parms &plist, int *ncrosscom);
int get_distance_2D(Fullmol *bases, double **traj, int p1, int p2, double deltat, double bindrad, int *ncross, int **crosspart, int **crossint, int i1, int i2, int mu, int **cross_rxn, int *ihome, double iter, double Rmax, double &sep, double &R1,  int *ncrosscom);
void get_distancePBCCELLMEM(int *numpartners, int **Speclist, int **myrxn, double mydensity, Protein *wholep, Complex *ind_com, double *kr, Fullmol *bases, double deltat, double *bindrad, int *ncross, int **crosspart, double **prob, int **crossint, int **cross_rxn, Parms &plist);
/*If particle has partners in reaction zone, evaluate whether a reaction will occur, and which one*/
int choose_one_reaction(double rnum, int p1, int *ncross, int **crosspart, double **prob, int &ci1, int &ci2, int **cross_rxn, int **cross_int, double irandmax);
int choose_one_reaction_loop(int p1, int *ncross, int **crosspart, double **prob, int &ci1, int &ci2, int **cross_rxn, int **cross_int, double irandmax);



/*rd input files, and checks on input files*/
void read_restart(ifstream &restart, int Ntotalmol, Fullmol *bases, Parms &plist, Complex *ind_com);
int read_restartCELL(ifstream &restart, int Ntotalmol, Fullmol *bases, Parms &plist, Complex *ind_com, int *Ncopy, Protein *wholep);
void read_parms(ifstream &parmfile, Parms &plist);
void read_parms_fits(ifstream &parmfile, Parms &plist, Fitparms &flist);
void read_parmsN(ifstream &parmfile, Parms &plist, std::vector<std::string> &infofilenames);
void read_reactions(ifstream &rxnfile, Parms &plist, int **Rlist, int *Ncoup, int **mycoupled, int *Nmyrxn, int **myrxn, double *bindrad, double *kr, int *cntrxn, int *freelist,int *bndlist, int *zlist, int *rxtype, double *Kd);
void update_rates(Protein *wholep, Parms &plist, int **Rlist,  double *bindrad, double *kr, int *rxtype, int clrxn, double *Kd, int *phome);
void check_reactions(Parms plist, Fullmol *bases, int *numpartners, int **Speclist, int *Nmyrxn, int **Rlist, int **myrxn, int Nspecies);
void set_status(ifstream &startfile, Protein *wholep, Fullmol *bases, Parms &plist, int *Ncopy);
void read_network(Parms &plist, int *numpartner, int **Speclist, ifstream &infile);
void read_protlist(int nwhole, Protein *wholep, int Niface, int *p_home, ifstream &protfile, int *ihome);
void renumber_list(int &nfree, int *freelist);
void read_coords(ifstream &crdfile, int Nprotypes, Protein *wholep, Fullmol* bases, Complex *ind_com, int *Ncopy, string *names);
void read_inputs(ifstream &rxnfile, Fullmol *bases, int *Ncopy, int *p_home, Protein *wholep, int *numpartner, int **Speclist, Parms &plist, int **Rlist, int *Ncoup, int **mycoupled, int *Nmyrxn, int **myrxn, double *bindrad, double *kr, int *cntrxn, int *freelist, int *bndlist, int *zlist, int *rxtype, int *ihome, std::vector<std::string> &infofilenames, double ***coordcont, int &howmanylipids);
void par_read_inputs(ifstream &rxnfile, Fullmol *bases, int *Ncopy, int *p_home, Protein *wholep, int *numpartner, int **Speclist, Parms &plist, int **Rlist, int *Ncoup, int **mycoupled, int *Nmyrxn, int **myrxn, double *bindrad, double *kr, int *cntrxn, int *freelist, int *bndlist, int *zlist, int *rxtype, int *ihome, std::vector<std::string> &infofilenames, double ***coordcont, int &howmanylipids);
void input_arrays(Fullmol *bases, int *Ncopy, int *p_home, Protein *wholep, int *numpartner, int **Speclist, Parms &plist, int **Rlist, int *Ncoup, int **mycoupled, int *Nmyrxn, int **myrxn, double *bindrad, double *kr, int *cntrxn, int *freelist, int *bndlist, int *zlist, int *rxtype, int *ihome);
double set_Rmaxlimit(Parms plist, Protein *wholep, double *bindrad, int *i_home, int *p_home, int **Rlist, int *rxtype, double ***coordcont, double &Rmax_radius);
double set_Rmaxlimit3_3(Parms plist, Protein *wholep, double *bindrad, int *i_home, int *p_home, int **Rlist, int *rxtype, double ***coordcont, double &Rmax_radius);


/*Enforce reflecting boundary at box perimeter*/
void reflect_traj_bound(int p1,  Fullmol *bases, Complex *ind_com, double xboxl, double yboxl, double zboxl, double **traj);
void reflect_traj_complex_rad_rot(int p1,  Fullmol *bases, Complex *ind_com, Parms plist, double *traj, double *trajR, double *M);
void reflect_traj_complex_rad_rot_nocheck(int p1,  Fullmol *bases, Complex *ind_com, double xboxl, double yboxl, double zboxl, double *traj, double *M);
void reflect_traj_check_span(double xtot, double ytot, double ztot, int p1, Fullmol *bases, Complex *ind_com, Parms plist, double *traj, double *trajR, double *M);
void reflect_traj_complex_rad(int p1,  Fullmol *bases, Complex *ind_com, double xboxl, double yboxl, double zboxl, double **traj);
void reflect_traj_complex_radCELL(int p1,  Fullmol *bases, Complex *ind_com, double xboxl, double yboxl, double zboxl, double **traj);
void reflect_bound(int p1,  Fullmol *bases, Complex *ind_com, double xboxl, double yboxl, double zboxl);
void reflect_complex_rad(int p1,  Fullmol *bases, Complex *ind_com, double xboxl, double yboxl, double zboxl);
void reflect_complex_rad_rot(int p1,  Fullmol *bases, Complex *ind_com, Parms plist);
void reflect_complex_rad_rotCELL(int p1,  Fullmol *bases, Complex *ind_com, double xboxl, double yboxl, double zboxl);
void reflect_traj_complex_rad_rotCELL(int p1,  Fullmol *bases, Complex *ind_com, double xboxl, double yboxl, double zboxl, double **traj, double *M);

/*Associate, and update all vectors*/
void perform_association(int p1, Fullmol *bases, Complex *ind_com, Parms &plist, int **crosspart, int **cross_rxn, int **crossint, int **Rlist, int *i_home, double **probvec, int *Ncoup, int **mycoupled, int *p_home, int ci1, int ci2, int it, int *movestat, int *ncross, int *ncrosscom);
void perform_association_sigma(int p1, Fullmol *bases, Complex *ind_com, Parms &plist, int **crosspart, int **cross_rxn, int **crossint, int **Rlist, int *i_home, double **probvec, int *Ncoup, int **mycoupled, int *p_home, int ci1, int ci2, int it, int *movestat, int *ncross, int bindrad, int *ncrosscom);
void perform_association_sigma_comPBCCELL(int p1, Fullmol *bases, Complex *ind_com, Parms &plist, int **crosspart, int **cross_rxn, int **crossint, int **Rlist, int *i_home, double **probvec, int *Ncoup, int **mycoupled, int *p_home, int ci1, int ci2, int it, int *movestat, int *ncross, int bindrad, int *ncrosscom, double **traj,int &newlegcount, int &newap2count, double **trajR);
void perform_association_sigma_com(int p1, Fullmol *bases, Complex *ind_com, Parms &plist, int **crosspart, int **cross_rxn, int **crossint, int **Rlist, int *i_home, double **probvec, int *Ncoup, int **mycoupled, int *p_home, int ci1, int ci2, int it, int *movestat, int *ncross, int bindrad, int *ncrosscom, double **traj,double **trajR, int *nBoundPairs);
void perform_association_sigma_com_crowd(int p1, Fullmol *bases, Complex *ind_com, Parms &plist, int **crosspart, int **cross_rxn, int **crossint, int **Rlist, int *i_home, double **probvec, int *Ncoup, int **mycoupled, int *p_home, int ci1, int ci2, int it, int *movestat, int *ncross, int bindrad, int *ncrosscom, double **traj,double **trajR, int *nBoundPairs);
void perform_association_sigma_com_crowd_stay(int p1, Fullmol *bases, Complex *ind_com, Parms &plist, int **crosspart, int **cross_rxn, int **crossint, int **Rlist, int *i_home, double **probvec, int *Ncoup, int **mycoupled, int *p_home, int ci1, int ci2, int it, int *movestat, int *ncross, int bindrad, int *ncrosscom, double **traj,double **trajR, int *nBoundPairs);
void perform_association_sigma_com_crowd_stayPBC(int p1, Fullmol *bases, Complex *ind_com, Parms &plist, int **crosspart, int **cross_rxn, int **crossint, int **Rlist, int *i_home, double **probvec, int *Ncoup, int **mycoupled, int *p_home, int ci1, int ci2, int it, int *movestat, int *ncross, int bindrad, int *ncrosscom, double **traj,double **trajR, int *nBoundPairs);
void perform_association_sigma_comSPHERES(int p1, Fullmol *bases, Complex *ind_com, Parms &plist, int **crosspart, int **cross_rxn, int **crossint, int **Rlist, int *i_home, double **probvec, int *Ncoup, int **mycoupled, int *p_home, int ci1, int ci2, int it, int *movestat, int *ncross, int bindrad, int *ncrosscom, double **traj);

/*Associate interfaces*/
void associate_norot(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist);
void associate_int_zfirst(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist);
void associate_int_zfirstPBC(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist);
int associate_freeleg(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom);
void associate_translate_measure(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist);
void associate_freelegPBC(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, int *ncrosscom);
void associate_translate_measurePBC(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, int *ncrosscom);
int associate_freelegPBCCELL(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, int *ncrosscom);
void associate_translate_measurePBCCELL(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, int *ncrosscom);
void associate_translate_measurePBCCELLMEM(int p1, int p2, int mu, int i1, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, Protein *wholep);
void associate_translate_measurePBCCELLMF(int p1, int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist,Protein *wholep, int *ncrosscom);
int associate_zlegs(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, int **inst_pro, int &mu_ret, int *knee_iface, int **myrxn, int *Nmyrxn);
void associate_zsigmaPBC(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom );
int associate_zsigmaPBCCELL(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom );
int associate_zsigma(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom );
int associate_crowd(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom );
int associate_crowd_stay(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom );
int associate_lipid_sigmaPBCCELL(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom );
int associate_lipid_sigma(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom );
int associate_lipid_sigmaPBCCELLAP2(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom);
int associate_lipid_sigmaAP2(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom);
int associate_ap2_clath_sigmaPBCCELL(int p1, int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom);
int associate_ap2_clath_sigma(int p1, int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom);
int associate_2dclat_3dclat_lipid_sigmaPBCCELL(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom);
int associate_2dclat_3dclat_lipid_sigma(int p1,int p2, int mu, int i1, int i2, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double bindrad, int *ncrosscom);


/*Update properties following association*/
void update_complex_propsPBC(int c1, Parms plist, Complex *ind_com, Fullmol *bases);
void update_complex_propsPBCCELL(int c1, Parms plist, Complex *ind_com, Fullmol *bases);
void update_complex_props(int c1, Parms plist, Complex *ind_com, Fullmol *bases);
void update_NofEach(int c1, Parms plist, Complex *ind_com, Fullmol *bases);
void update_bound_proteinsPBC(int p1, int p2, int i1, int i2, Fullmol *bases,  Complex *ind_com, Parms &plist, int *ncrosscom, int prod, int iind, int iind2);
void update_bound_proteinsPBCCELL(int p1, int p2, int i1, int i2, Fullmol *bases,  Complex *ind_com, Parms &plist, int *ncrosscom, int prod, int iind, int iind2);
void update_bound_proteins(int p1, int p2, int i1, int i2, Fullmol *bases,  Complex *ind_com, Parms &plist, int *ncrosscom, int prod, int iind, int iind2);
void update_crowded_proteins(int p1, int p2, int i1, int i2, Fullmol *bases,  Complex *ind_com, Parms &plist, int *ncrosscom, int prod, int iind, int iind2);
void update_Nboundpairs(int ptype1, int ptype2, int chg,  Parms plist, int *nBoundPairs);



/*Dissociate interfaces*/
int break_complex_trans(int p1, int mu, int kind, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double *bindrad,  int p2, int i1, int i2, int *p_home, int **myrxn);
int break_complexPBC(int p1, int mu, int kind, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double *bindrad,  int p2, int i1, int i2, int *p_home, int **myrxn );
int break_sigmaPBC(int p1, int mu, int kind, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double *bindrad,  int p2, int i1, int i2, int *p_home, int **myrxn);
int break_sigmaPBCCELL(int p1, int mu, int kind, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double *bindrad,  int p2, int i1, int i2, int *p_home, int **myrxn);
int break_sigma(int p1, int mu, int kind, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double *bindrad,  int p2, int i1, int i2, int *p_home, int **myrxn);
int break_complexPBCCELL(int p1, int mu, int kind, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double *bindrad,  int p2, int i1, int i2, int *p_home, int **myrxn );
int break_complexPBCCELLMEM(int p1, int mu, int kind, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double *bindrad, int p2, int i1, int i2, int *p_home, int **myrxn,Protein *wholep);
int determine_which_complex_double(int p1, int p2, Complex *ind_com, int c1, int c2, Fullmol *bases, int **myrxn, int **Rlist, int *ihome, int *p_home, Parms &plist);
int break_complexPBCCELLMF(int p1, int mu, int kind, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double *bindrad, int p2, int i1, int i2, int *p_home, int **myrxn,Protein *wholep);
int break_complex_zleg_clath(int p1, int mu, int kind, Fullmol *bases, int **Rlist, int *ihome, Complex *ind_com, Parms &plist, double *bindrad, int **inst_pro, int p2, int clathstop, int **myrxn, int *Ncoup, int **mycoup, int *p_home);



/*Special calls for clathrin trimer*/
int get_kneeind_canc(Fullmol *bases, int piv2, int orig2, int *knee_iface, int *ihome);
int get_kneeind2(Fullmol *bases, int piv2, int orig2, int *knee_iface, int *ihome);
int measure_overlap(int c1, int c2, Complex *ind_com, Fullmol *bases, int pclath);

/*For reflecting BC, check for span*/
int check_box_span(int p1, int p2, Fullmol *bases, Complex *ind_com, Parms &plist);

/*preventing large movement to bind within the same complex*/
int same_complex_test(int p1, int p2, Fullmol *bases, int **crossint, int *i_home, double *bindrad, int rxn1, int ci1, int ci2);


/*Calling coupled reactions*/
void ap2_coupledrxn_add(int mu, int *Nmycoup, int **mycoupled, Fullmol *bases, int **Rlist , int *i_home, int *p_home, int p1, int p2);
void ap2_coupledrxn_sub(int mu, int *Nmycoup, int **mycoupled, Fullmol *bases, int **Rlist , int *i_home, int *p_home, int p1, int p2);
void update_coupledrxn_add(int mu, int *Nmycoup, int **mycoupled, Fullmol *bases, int **Rlist, int *i_home, int *p_home, int p1, int p2 );


/*Update particle positions and also COM for some*/
void move_protein_pos(int p1,  Fullmol *bases, Complex *ind_com, double **traj, int *movestat);
void move_protein_posPBC(int p1,  Fullmol *bases, Complex *ind_com, double **traj, int *movestat, Parms &plist);
void move_zcrds(Fullmol *bases, int p1, int p2, double delz1, double delz2, Complex *ind_com);
void move_zcrdsPBC(Fullmol *bases, int p1, int p2, double delz1, double delz2, Complex *ind_com, Parms &plist);
void move_compro_coords(int c1, Complex *ind_com, Fullmol *bases, double *chg1);
void move_pro_coordsPBC(int c1, Complex *ind_com, Fullmol *bases, double *chg1, Parms &plist);
void move_pro_coordsPBCCELL(int c1, Complex *ind_com, Fullmol *bases, double *chg1, Parms &plist);
void move_pro_coords(int c1, Complex *ind_com, Fullmol *bases, double *chg1);

/*rotation calls for proteins*/
void rotate_and_translate(int p1, int c1, Complex *ind_com, Fullmol *bases, double *M, double *dtrans);
void rotate_and_translate_int(int p1, int c1, Complex *ind_com, Fullmol *bases, double *M, double *dtrans, int iind);
void move_rot_proteins(int p1,  Fullmol *bases, Complex *ind_com, double *traj, int *movestat, double *trajR, double *M);
void move_rot_proteinsPBC(int p1,  Fullmol *bases, Complex *ind_com, double *traj, int *movestat, double *trajR, double *M, Parms plist);
void move_rot_proteinsPBC_substep(int p1,  Fullmol *bases, Complex *ind_com, double *traj, int *movestat, double *trajR, double *M, Parms plist);
void rot_euler_1p(int p1, Fullmol *bases, double deltat, double *M);
void rotate_only(int p1, int c1, Complex *ind_com, Fullmol *bases, double *M);
void rotate_onlyPBC(int p1, int c1, Complex *ind_com, Fullmol *bases, double *M, Parms plist);
void translate(double *v, double *d, double *v2);
void translate_int(int p1, int c1, Complex *ind_com, Fullmol *bases, double *dtrans);
void translate_intPBC(int p1, int c1, Complex *ind_com, Fullmol *bases, double *dtrans, Parms &plist);
void rotate_and_translate_intPBC(int p1, int c1, Complex *ind_com, Fullmol *bases, double *M, double *dtrans, int iind, Parms &plist);
void move_rot_proteinsPBCCELL(int p1,  Fullmol *bases, Complex *ind_com, double **traj, int *movestat, double **trajR, double *M, Parms plist);
void rotate_onlyPBCCELL(int p1, int c1, Complex *ind_com, Fullmol *bases, double *M, Parms plist);
void translate_intPBCCELL(int p1, int c1, Complex *ind_com, Fullmol *bases, double *dtrans, Parms &plist);
void rotate_and_translate_intPBCCELL(int p1, int c1, Complex *ind_com, Fullmol *bases, double *M, double *dtrans, int iind, Parms &plist);

/*move particles and avoid others, but doesn't guarantee that there won't be overlap at end*/
void move_and_exclude_other(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad);
void move_and_exclude_seq(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad);
void move_and_exclude_rot(double deltat, int p1, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, int **cross_int, int **cross_rxn, double **traj, double **prob, Parms plist, int *movestat, int *ihome, double *bindrad, double *M, double **trajR);



/*After association or dissociation, change D and radius*/
void update_diffusion(int c1, Complex *ind_com, Fullmol *bases);
void update_diffusionMF(int c1, Complex *ind_com, Fullmol *bases,Protein *wholep);
void update_radius(int c1, Complex *ind_com, Fullmol *bases);
void update_radiusPBC(int c1, Complex *ind_com, Fullmol *bases, Parms plist);
void update_radius_rot(int c1, Complex *ind_com,  Fullmol *bases);
void update_radiusPBCCELL(int c1, Complex *ind_com, Fullmol *bases, Parms plist);

/*updating diffusion using Einstein-stokes scale factors and radius, rather than
  D read in from file.*/
double rot_prefactor(double T, double nu, double scale, double a, double Drot);
double trans_prefactor(double T, double nu, double scale, double a, double Dtot);
void update_trans_diffusion(int c1, Complex *ind_com, Fullmol *bases, double pretrans);
void update_rot_diffusion(int c1, Complex *ind_com, Fullmol *bases);
void update_rot_diffusion_a3(int c1, Complex *ind_com, Fullmol *bases, double prerot) ;


/*write out coordinates, etc.*/
void write_status(ofstream &outfile, Protein *wholep, Fullmol *bases, Parms &plist, int *Ncopy, int iter);
void write_statusCELL(ofstream &outfile, ofstream &outfile2, Protein *wholep, Fullmol *bases, Parms &plist, int *Ncopy, int iter, Complex *ind_com, double deltat);
void write_complex_size(int Nc, Complex *ind_com, Fullmol *bases, ofstream &outfile ,int it);
void write_complex(ofstream &outfile, Parms &plist, Complex *ind_com, double time);
void write_crds(Fullmol *bases, int p1);
void write_crds_complex(int c1, Complex *ind_com, Fullmol *bases);
void write_separation(Fullmol *bases, int p1, int p2);
void write_inst(ofstream &outfile, Parms &plist, int **inst_pro, Protein *wholep, int it);
void write_protein(ofstream &outfile, Parms &plist, Fullmol *bases, int *Ncopy, int it);
void write_protein_iface(ofstream &outfile, Parms &plist, Fullmol *bases, int *Ncopy, int it, Protein *wholep);
void write_protein_iface_short(ofstream &outfile, Parms &plist, Fullmol *bases, int *Ncopy, int it, Protein *wholep, string *names);
void write_dcd(ofstream &outfile, Parms &plist, Fullmol *bases, int *Ncopy, int it, Protein *wholep, std::vector<std::string> &infofilenames);
void write_protein_iface_short_stretch(ofstream &outfile, Parms &plist, Fullmol *bases, int *Ncopy, int it, Protein *wholep, string *names, double stretch);
void copy_crds(int Nmol, Fullmol *bases, double *crd);
void copy_iface_crds(int Nmol, Fullmol *bases, double *crd, Protein *wholep);
void write_complex_cout(int Nc, Complex *ind_com, Fullmol *bases);
void write_timestat(ofstream &outfile, Protein *wholep, Fullmol *bases, Parms &plist, int *Ncopy, int iter, Complex *ind_com, double deltat,int Nprotypes);
void write_timestat2(ofstream &outfile, ofstream &molectypesfile, ofstream &timestatfiletext, Protein *wholep, Fullmol *bases, Parms &plist, int *Ncopy, int iter, Complex *ind_com, double deltat,int Nprotypes, int **molectypes, int &currentnumberofmolectypes);
void write_restart(ofstream &outfile, Protein *wholep, Fullmol *bases, Parms &plist, int *Ncopy, int iter, Complex *ind_com, double deltat);
void gen_psf_system(Parms &plist, Protein *wholep, int *Ncopy,std::vector<std::string> &infofilenames);

double calc_complex_hist(int Nc, Complex *ind_com, Fullmol *bases, ofstream &outfile, int it, Parms plist, std::vector<string> &infofilenames, int *Ncopy);
void print_dimers(int Nc, Complex *ind_com, Fullmol *bases, ofstream &outfile, int it, Parms plist, std::vector<string> &infofilenames, int *Ncopy);
void init_print_dimers( ofstream &outfile, int it, Parms plist, std::vector<string> &infofilenames);
void init_NboundPairs(int *nBoundPairs, std::vector<int> &proPairlist, ofstream &outfile,  Parms plist, std::vector<string> &infofilenames, Protein *wholep);
void write_NboundPairs(int *nBoundPairs, std::vector<int> &proPairlist,  ofstream &outfile, int it, Parms plist);

void write_complex_components(int Nc, Complex *ind_com, Fullmol *bases, ofstream &outfile, int it, Parms plist, std::vector<string> &infofilenames);


/*calculate gofrs for nonPBC boxes*/
void calc_gr(Parms &plist, double delr, int nbins, Fullmol *bases, int *Ncopy, double **gr, double **Vofr);
void calc_gr_center(Parms plist, double delr, int nbins, Fullmol *bases,double *gr, double bindrad);
void calc_gr_self(Parms plist, double delr, int nbins, Fullmol *bases,double *gr, double bindrad);
void calc_gr_self_nnorm(Parms plist, double delr, int nbins, Fullmol *bases,double *gr, double bindrad);
void init_gr_zero(double *gr, int nbins);
void gr_volume_norm(double *gr, int nbins, double delr, double rho0, int Nrep, double bindrad, int flagdim);

void write_gofr2(int n1, int n2, Parms &plist, int nbins, double **gr, double delr, ofstream &outfile);
void write_gofr(int nbins, double *gr, double delr, ofstream &outfile);
void get_volumes(int nbins, Fullmol *bases, double **Vofr, Parms plist,  double delr, int N);

/*sample positions with PBC*/
void generate_initial_crds_ABavoid(Parms plist, Fullmol *bases,int *Ncopy, Complex *ind_com, double *bindrad, int zeroB);
void generate_initial_crds_AB_allavoid(Parms plist, Fullmol *bases,int *Ncopy, Complex *ind_com, double *bindrad);
void gen_copy_coords(int Nprotypes, Protein *wholep, Fullmol* bases, Complex *ind_com, int *Ncopy);
void generate_initial_crds_AB(Parms plist, Fullmol *bases,int *Ncopy, Complex *ind_com, double *bindrad);
/*void generate_initial_crds(Parms plist, Fullmol *bases, int *Ncopy, Complex *ind_com, double *bindrad, int Nprotypes, int Nifaces, Protein *wholep, int **Rlist,int *rxtype, int *p_home);*/
void generate_initial_crdsPBC(Parms plist, Fullmol *bases, int *Ncopy, Complex *ind_com, double *bindrad, Protein *wholep, int **Rlist,int *rxtype, int *p_home, int howmanylipids, double ***coordcont);
void generate_initial_crdsPBCCELL(Parms plist, Fullmol *bases, int *Ncopy, Complex *ind_com, double *bindrad, Protein *wholep, int **Rlist,int *rxtype, int *p_home, int howmanylipids, double ***coordcont);
void generate_initial_crdsNoPBC(Parms plist, Fullmol *bases, int *Ncopy, Complex *ind_com, double *bindrad, Protein *wholep, int **Rlist,int *rxtype, int *p_home, int howmanylipids, double ***coordcont);
void generate_initial_crds_legs(Parms plist, Fullmol *bases, int *Ncopy, Complex *ind_com, double *bindrad, Protein *wholep, int **Rlist,int *rxtype, int *p_home, int howmanylipids, double ***coordcont);

/*debug calls*/
void check_movement(int Nmol, Fullmol *bases, double *crd);
void check_com_to_pro(int Nmol, int Ncomplex, Fullmol *bases, Complex *ind_com);
void check_pro_to_com(int Nmol, Fullmol *bases, Complex *ind_com);
void check_bndlist(int Ncomplex, Complex *ind_com, Fullmol *bases);
void check_membrane(int *Ncopy, Fullmol *bases);
void check_rot_move(int Nmol, Fullmol *bases, Complex *ind_com);

/*calls for breaking that are not working*/
void break_and_avoid(int pB, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, double **traj, Parms plist, double bindrad);
void break_and_avoid_both(int pB, int pA, Fullmol *bases, Complex *ind_com, int *ncross, int **cross_part, double **traj, Parms plist, double bindrad);
void get_overlap_list(int p1, Fullmol *bases, int *ncross, int **crosspart, double bindrad, int *binlist, int MAXPERBIN,int maxnbor, int *npb, int *nbor, int *nborrev, Parms plist, double maxsep2);


/*old methods to determine collisions, and pick partners*/
void integrate_path(Fullmol *bases, double **traj, int p1, int p2, double deltat, double bindrad, int *ncross, double **treturn, int **crosspart, int **crossint, int i1, int i2, int mu, int **cross_rxn, int *ihome, double iter);
void test_path(int j, Fullmol *bases, int wprot, Protein *wholep, int i, int *numpartners, int **Speclist, int **myrxn, double *bindrad, double deltat, double **traj, int *ncross, int *i_home, double **treturn, int **crosspart, int **crossint, int **cross_rxn);
void move_prod_pos(double tremain, int p1,  int p2, int i1, int i2, Fullmol *bases, Complex *ind_com);
int choose_large_reaction(double rnum, int p1, int *ncross, int **crosspart, double **prob, int &ci1, int &ci2, int **cross_rxn, int **cross_int);
int choose_one_partner(int p1, int *ncross, int **crosspart, double **treturn, int &ci1, int &ci2);
int choose_one_partner_ofm(int p1, int *ncross, int **crosspart, double **treturn, int &ci1, int &ci2, int **cross_rxn, int **cross_int);





#endif
