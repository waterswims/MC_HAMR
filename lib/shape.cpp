#include "../includes/shape.hpp"

#ifdef __INTEL_COMPILER
#include "../includes/mklrand.hpp"
#define DRANDTYPE mklrand::mkl_drand
#else
#include "../includes/stdrand.hpp"
#define DRANDTYPE stdrand::std_d_unirand
#endif

#include <cmath>
#include <iostream>

extern DRANDTYPE st_rand_double;

particle::shape::weibull::weibull(shape_type& other)
{
    r0 = other.get_r0();
    beta = other.get_beta();
    a[0] = other.get_a();
    a[1] = other.get_b();
    a[2] = other.get_c();
}

particle::shape::weibull::weibull(double rin, double bin)
{
    beta = bin;
    r0 = 1 / tgamma(1. + 1. / beta);
    a[0] = rin;
    a[1] = rin;
    a[2] = rin;
}

particle::shape::weibull::weibull(double betain,
    double ain,
    double bin,
    double cin)
{
    beta = betain;
    r0 = 1 / tgamma(1. + 1. / beta);
    a[0] = ain;
    a[1] = bin;
    a[2] = cin;
}

bool particle::shape::weibull::check(std::vector<int> Is, int l_size)
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

particle::shape::weibull& particle::shape::weibull::operator=(shape_type& other)
{
    r0 = other.get_r0();
    beta = other.get_beta();
    a[0] = other.get_a();
    a[1] = other.get_b();
    a[2] = other.get_c();
    return *this;
}
