#ifndef _FUNC
#define _FUNC

#include <vector>
#include <fstream>
#include <cstring>

#define XTENSOR_USE_XSIMD
#include <xtensor/xfixed.hpp>

double mean(std::vector<double> &oY);

double std_dev(std::vector<double> &x);

double norm(std::vector<double> vals);

double sum(std::vector<double> &oY);

int sum(std::vector<int> &oY);

void AtoLn(double amean, double asd, double &lmean, double &lsd);

int mod(int a, int b);

xt::xtensorf<double, xt::xshape<4>> c_prod(
    const xt::xtensorf<double, xt::xshape<4>> &s1,
    const xt::xtensorf<double, xt::xshape<4>> &s2);

#endif
