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

double solid_angle(const std::vector<double> &s1,
                const std::vector<double> &s2,
                const std::vector<double> &s3,
                std::vector<double> &buff)
{
    if(s1[0] == 0 && s1[1] == 0 && s1[2] == 0) {return 0;}
    else if(s2[0] == 0 && s2[1] == 0 && s2[2] == 0) {return 0;}
    else if(s3[0] == 0 && s3[1] == 0 && s3[2] == 0) {return 0;}

    double s1s2 = 0, s1s3 = 0, s2s3 = 0;
    #pragma omp simd
    for(unsigned int i = 0; i < 4; i++)
    {
        s1s2 += s1[i] * s2[i];
        s1s3 += s1[i] * s3[i];
        s2s3 += s2[i] * s3[i];
    }
    double rho = 2 * (1 + s1s2) * (1 + s1s3) * (1 + s2s3);
    double dotsum = 1 + s1s2 + s1s3 + s2s3;
    buff[0] = s1[1]*s3[2] - s1[2]*s3[1];
    buff[1] = s1[2]*s3[0] - s1[0]*s3[2];
    buff[2] = s1[0]*s3[1] - s1[1]*s3[0];
    double crosssum = 0;
    #pragma omp simd
    for(unsigned int i = 0; i < 4; i++)
    {
        crosssum += buff[i] * s2[i];
    }
    double cossum = dotsum / rho;
    double sinsum = crosssum / rho;
    double ang = atan2(sinsum, cossum);

    return ang;
}

std::valarray<double> c_prod(const std::valarray<double> &s1,
    const std::valarray<double> &s2)
{
    std::valarray<double> out(s1);
    out[0] = s1[1]*s2[2] - s1[2]*s2[1];
    out[1] = s1[2]*s2[0] - s1[0]*s2[2];
    out[2] = s1[0]*s2[1] - s1[1]*s2[0];

    return out;
}
