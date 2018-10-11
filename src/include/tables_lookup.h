#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>

using namespace std;

double look_up(double rnum, double delp, double *table);
double look_up_dr(double rnum, double *delp, double **table, double sep, double delR);
void create_table(double maxprob, double *table, double MAXP, double kact, double kdiff, double Dtot, double bindrad);
void create_table_dt(double maxprob, double *table, double MAXP, double kact, double kdiff, double Dtot, double bindrad, double R1, double alpha, double &delp);
void create_2D_table(int MAXR, double Dtot, double bindrad, double **table, double *delp, double delR, double kact, double kdiff, int MAXP, double deltat);
