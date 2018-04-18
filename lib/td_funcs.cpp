#include "../includes/td_funcs.hpp"
#include "../includes/functions.hpp"

#include <cmath>
#include <iostream>

double particle::funcs::calc_E(particle::field::field_type& lattice,
    std::valarray<double> H)
{
    int Nspins = lattice.get_size();
    std::valarray<double> J_sum(lattice.access(0).size()),
        H_sum(lattice.access(0).size()), D_sum(lattice.access(0).size());

    J_sum = 0;
    H_sum = 0;
    D_sum = 0;

    for(int i = 0; i < Nspins; i++)
    {
        for(int j = 0; j < lattice.get_neigh(i).size(); j++)
        {
            J_sum += lattice.get_J(i, j) *
                lattice.access(lattice.get_neigh(i)[j]) *
                lattice.access(i);

            D_sum += lattice.get_D_vec(i, j) *
                c_prod(lattice.access(lattice.get_neigh(i)[j]),
                lattice.access(i));
        }
    }
    H_sum = calc_M(lattice);

    double E = (-H*H_sum - 0.5*J_sum - 0.5*D_sum).sum();

    return E;
}

double particle::funcs::calc_dE(particle::field::field_type& lattice,
    int position, std::valarray<double> H)
{
    lattice.gen_rand();
    std::valarray<double> diff = lattice.access(position) - lattice.get_rand();
    std::valarray<double> J_sum(diff.size()), D_sum(diff.size());
    J_sum = 0;
    D_sum = 0;

    for(int j = 0; j < lattice.get_neigh(position).size(); j++)
    {
        J_sum += lattice.get_J(position, j) *
            lattice.access(lattice.get_neigh(position)[j]);
        D_sum += c_prod(lattice.access(lattice.get_neigh(position)[j]),
            lattice.get_D_vec(position, j));
    }

    double dE = ((H + J_sum + D_sum) * diff).sum();

    return dE;
}

std::valarray<double> particle::funcs::calc_M(
    particle::field::field_type& lattice)
{
    int Nspins = lattice.get_size();
    std::valarray<double> s(lattice.access(0).size());
    s = 0;

    for(int i = 0; i < Nspins; i++)
    {
        s += lattice.access(i);
    }

    return s;
}

std::valarray<double> particle::funcs::calc_subM(
    particle::field::field_type& lattice,
    int subnumber)
{
    int Nspins = lattice.get_size();
    std::valarray<double> s(lattice.access(0).size());
    s = 0;

    for(int i = 0; i < Nspins; i++)
    {
        int possum = lattice.get_loc(i).sum();
        if (possum%2 == subnumber)
        {
            s += lattice.access(i);
        }
    }

    return s;
}

std::vector<double> particle::funcs::calc_TC(
    particle::field::field_type& lattice)
{
    int Nspins = lattice.get_size();
    std::vector<double> out;
    std::valarray<int> loc_diff;
    int left, right, forward, backward;

    for(int i = 0; i < Nspins; i++)
    {
        if(out.size() < lattice.get_loc(i)[2]+1)
        {
            out.resize(out.size());
        }
        for(int j = 0; j < Nspins; j++)
        {
            loc_diff = lattice.get_loc(j) - lattice.get_loc(i);
            if(loc_diff[2] == 0)
            {
                if (loc_diff[0] == 1 && loc_diff[1] == 0)
                {
                    left = j;
                }
                else if (loc_diff[0] == -1 && loc_diff[1] == 0)
                {
                    right = j;
                }
                else if (loc_diff[0] == 0 && loc_diff[1] == 1)
                {
                    forward = j;
                }
                else if (loc_diff[0] == 0 && loc_diff[1] == -1)
                {
                    backward = j;
                }
            }
        }
        out[lattice.get_loc(i)[2]] += solid_angle(lattice.access(i),
            lattice.access(right), lattice.access(forward));
        out[lattice.get_loc(i)[2]] += solid_angle(lattice.access(i),
            lattice.access(left), lattice.access(backward));
    }
    for(int i = 0; i < out.size(); i++)
    {
        out[i] /= 2 * M_PI;
    }

    return out;
}

double particle::funcs::solid_angle(const std::valarray<double> &s1,
                const std::valarray<double> &s2,
                const std::valarray<double> &s3)
{
    double s1s2 = (s1 * s2).sum();
    double s1s3 = (s1 * s3).sum();
    double s2s3 = (s2 * s3).sum();

    double rho = 2 * (1 + s1s2) * (1 + s1s3) * (1 + s2s3);
    double dotsum = 1 + s1s2 + s1s3 + s2s3;

    std::valarray<double> s1xs3 = c_prod(s1, s3);
    double crosssum = (s1xs3 * s2).sum();

    double ang = atan2(crosssum / rho, dotsum / rho);

    return ang;
}
