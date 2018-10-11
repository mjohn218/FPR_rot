#ifndef __CLASS_DEFS_H
#define __CLASS_DEFS_H

#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include <string>
#include <sys/time.h>

#define MAXIFACE 8
#define MAXPRTNER 10
#define MAXCOMPLEX 700
#define MAXRXN 200
#define MAXOVERLAP 100
#define MAXNUMCOMPLEXTYPES 1000
#define MAXNUMCELL 500000
#define KABS 1e15 //If K>KABS (nm^3/us) use absorbing BC solutions of GFs

using namespace std;

class Mbind
{
public:
  int p1start;
  int p2start;
  int ncol;
  int nrow;
};
class Protein
{
public:
  int ninterface;
  int valiface[MAXIFACE];
  int npropart;
  int propart[MAXPRTNER];
  double Dx;
  double Dy;
  double Dz;
  double Drx;
  double Dry;
  double Drz;
  double radR;
  //  double radx;
  //double rady;
  //double radz;
  int nint_write;
  int wrlist[MAXIFACE];
  //  vector<int> Npairs;
};
class Fullmol
{
public:
  double mytime;
  int mybin;
  int mybinind;
  int protype;
  int ninterface;
  int istatus[MAXIFACE];
  int npartner;
  int mycomplex;
  
  double xcom;
  double ycom;
  double zcom;
  double x[MAXIFACE];
  double y[MAXIFACE];
  double z[MAXIFACE];
  int nbnd;
  int nfree;
  int freelist[MAXIFACE];
  int bndlist[MAXIFACE];
  double mass;
  //double massx;
  //double massy;
  //double massz;
  int partner[MAXIFACE];
  double Dx;
  double Dy;
  double Dz;
  
  double radR;
  double Drx;
  double Dry;
  double Drz;
  int npropart;
  int propart[MAXPRTNER];
  //  double radx;
  //double rady;
  //double radz;
};
class Complex
{
public:
  int mysize;
  double Dx;
  double Dy;
  double Dz;
  int plist[MAXCOMPLEX];
  double xcom;
  double ycom;
  double zcom;
/*   double radxn; */
/*   double radyn; */
/*   double radzn; */
/*   double radxp; */
/*   double radyp; */
/*   double radzp; */
  double radR;
  double Drx;
  double Dry;
  double Drz;
  vector<int> NofEach;
 
};
class Parms
{
public:
  int Nprotypes;
  int Nifaces;
  double Nit;
  int restart;
  int statwrite;
  int configwrite;
  int grwrite;//This should be renamed restartWrite
  int Nrep;
  int grfreq;
  double deltat;
  double V;
  double X0total;
  int Nspecies;
  int Nrxn;
  double mass;
  double D;
  int nspec_complex;
  double maxsep2;
  int ntotalcomplex;
  double xboxl;
  double yboxl;
  double zboxl;
  int Ntotalmol;
  int Natom;
  int Natomwrite;
  double pretrans;
  double prerot;
  int pclath;
  int nloop;
};
class Fitparms
{
 public:
  double Dfit;
  double kfit;
  double sfit;

};

class ClusterPair
{

 public:
  int p1;
  int p2;
  int i1;
  int i2;
  int k1;
  int k2;
  int priority;
  int memtest;
  double bindrad;
  //int flagOverlap;

  ClusterPair(void); //Constructor
  ClusterPair(int setp1, int setp2); //Constructor
};

#endif
