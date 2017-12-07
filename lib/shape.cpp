#include "../includes/shape.hpp"
#include "../includes/mklrand.hpp"
#include <cmath>
#include <iostream>

extern mkl_drand st_rand_double;

weibull::weibull(shape_type& other)
{
    r0 = other.get_r0();
    beta = other.get_beta();
    a[0] = other.get_a();
    a[1] = other.get_b();
    a[2] = other.get_c();
}

weibull::weibull(double rin, double bin)
{
    beta = bin;
    r0 = 1 / tgamma(1. + 1. / beta);
    a[0] = rin;
    a[1] = rin;
    a[2] = rin;
}

weibull::weibull(double betain, double ain, double bin, double cin)
{
    beta = betain;
    r0 = 1 / tgamma(1. + 1. / beta);
    a[0] = ain;
    a[1] = bin;
    a[2] = cin;
}

bool weibull::check(std::vector<int> Is, int l_size)
{
    double centre = double(l_size - 1) / 2.;
    double dist2 = 0;
    for(int i = 0; i < Is.size(); i++)
    {
        dist2 += pow((Is[i]-centre)/(a[i]), 2);
    }
	double dist = pow(dist2, 0.5);
	double test = exp(-pow((dist/r0), beta));
	if(st_rand_double.gen() < test)
	{
		return true;
	}
	else
	{
		return false;
	}
}

weibull& weibull::operator=(shape_type& other)
{
    r0 = other.get_r0();
    beta = other.get_beta();
    a[0] = other.get_a();
    a[1] = other.get_b();
    a[2] = other.get_c();
    return *this;
}
