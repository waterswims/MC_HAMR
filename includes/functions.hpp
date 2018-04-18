#ifndef _FUNC
#define _FUNC

#include <vector>
#include <valarray>
#include <fstream>
#include <cstring>

double mean(std::vector<double> &oY);

double std_dev(std::vector<double> &x);

double norm(std::vector<double> vals);

double sum(std::vector<double> &oY);

int sum(std::vector<int> &oY);

void AtoLn(double amean, double asd, double &lmean, double &lsd);

int mod(int a, int b);

double solid_angle(const std::vector<double> &s1,
                const std::vector<double> &s2,
                const std::vector<double> &s3,
                std::vector<double> &buff);

std::valarray<double> c_prod(const std::valarray<double> &s1,
    const std::valarray<double> &s2);

#endif
