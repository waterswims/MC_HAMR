#ifndef _FUNC
#define _FUNC

#include <vector>
#include <fstream>
#include <cstring>

using namespace std;

double mean(vector<double> &oY);

double mean(double oY[]);

double std_dev(vector<double> &x);

double mag_sus(vector<double> magz, double T, double kb);

double heat_cap(vector<double> ener, double T, double kb);

double binders(vector<double> mag, int size);

int arr_len(double oY[]);

double norm(double val, int norm_typ);

double norm(vector<double> vals);

double sum(vector<double> &oY);

int sum(vector<int> &oY);

void printv(vector<int> &oX);

void AtoLn(double amean, double asd, double &lmean, double &lsd);

int mod(int a, int b);

#endif
