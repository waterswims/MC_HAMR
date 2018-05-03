#include "../includes/functions.hpp"
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <sstream>

double mean(std::vector<double> &oY){
    double s(0);
    for (int i=0; i < oY.size(); i++){
        s += oY[i];
    }
    return s / oY.size();
}

double std_dev(std::vector<double> &x){
    double av_x = mean(x);
    int size = x.size();
    double sum = 0;
    for(int i(0); i < size; i++){
        sum += pow(x[i] - av_x, 2);
    }
    double var = sum/(size - 1.);
    return pow(var, 0.5);
}

double norm(std::vector<double> vals)
{
    double s = 0;
    for(std::vector<double>::iterator it = vals.begin(); it != vals.end(); it++)
    {
        s += pow(*it, 2);
    }
    return pow(s, 0.5);
}

double sum(std::vector<double> &oY)
{
            double s = 0;
            for(std::vector<double>::iterator it = oY.begin(); it != oY.end(); it++)
            {
                        s += *it;
            }
            return s;
}

int sum(std::vector<int> &oY)
{
    int s = 0;
    for(std::vector<int>::iterator it = oY.begin(); it != oY.end(); it++)
    {
                s += *it;
    }
    return s;
}

void AtoLn(double amean, double asd, double &lmean, double &lsd)
{
    double av = pow(asd, 2);
    double am2 = pow(amean, 2);
    lmean = log(am2 / pow(av+am2, 0.5));
    lsd = pow(log(1 + av / am2), 0.5);
}

int mod(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}

xt::xtensorf<double, xt::xshape<4>> c_prod(
    const xt::xtensorf<double, xt::xshape<4>> &s1,
    const xt::xtensorf<double, xt::xshape<4>> &s2)
{
    xt::xtensorf<double, xt::xshape<4>> out;
    out[0] = s1[1]*s2[2] - s1[2]*s2[1];
    out[1] = s1[2]*s2[0] - s1[0]*s2[2];
    out[2] = s1[0]*s2[1] - s1[1]*s2[0];

    return out;
}
