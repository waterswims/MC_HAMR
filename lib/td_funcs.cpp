#include "../includes/td_funcs.hpp"
#include "../includes/functions.hpp"

#include <cmath>
#include <iostream>

#include <xtensor/xio.hpp>

xt::xtensorf<double, xt::xshape<4>> J_sum, H_sum, D_sum, diff;

double particle::funcs::calc_E(field::field_type& lattice,
    xt::xtensorf<double, xt::xshape<4>>& H)
{
    int Nspins = lattice.get_size();

    J_sum = {0, 0, 0, 0};
    H_sum = {0, 0, 0, 0};
    D_sum = {0, 0, 0, 0};

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

    double E = xt::sum(-H * H_sum - 0.5 * J_sum - 0.5 * D_sum)[0];

    return E;
}

double particle::funcs::calc_dE(particle::field::field_type& lattice,
    int position, xt::xtensorf<double, xt::xshape<4>>& H)
{
    lattice.gen_rand();
    diff = lattice.access(position) - lattice.get_rand();
    J_sum = {0, 0, 0, 0};
    D_sum = {0, 0, 0, 0};

    for(int j = 0; j < lattice.get_neigh(position).size(); j++)
    {
        J_sum += lattice.get_J(position, j) *
            lattice.access(lattice.get_neigh(position)[j]);
        D_sum += c_prod(lattice.access(lattice.get_neigh(position)[j]),
            lattice.get_D_vec(position, j));
    }

    double dE = xt::sum((H + J_sum + D_sum) * diff)[0];

    return dE;
}

xt::xtensorf<double, xt::xshape<4>> particle::funcs::calc_M(
    particle::field::field_type& lattice)
{
    int Nspins = lattice.get_size();
    xt::xtensorf<double, xt::xshape<4>> s;
    s = {0, 0, 0, 0};

    for(int i = 0; i < Nspins; i++)
    {
        s += lattice.access(i);
    }

    return s;
}

xt::xtensorf<double, xt::xshape<4>> particle::funcs::calc_subM(
    particle::field::field_type& lattice,
    int subnumber)
{
    int Nspins = lattice.get_size();
    xt::xtensorf<double, xt::xshape<4>> s;
    s = {0, 0, 0, 0};

    for(int i = 0; i < Nspins; i++)
    {
        int possum = xt::sum(lattice.get_loc(i))[0];
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
    xt::xtensorf<int, xt::xshape<4>> loc_diff;
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

double particle::funcs::solid_angle(const xt::xtensorf<double, xt::xshape<4>>& s1,
                const xt::xtensorf<double, xt::xshape<4>>& s2,
                const xt::xtensorf<double, xt::xshape<4>>& s3)
{
    double s1s2 = xt::sum(s1 * s2)[0];
    double s1s3 = xt::sum(s1 * s3)[0];
    double s2s3 = xt::sum(s2 * s3)[0];

    double rho = 2 * (1 + s1s2) * (1 + s1s3) * (1 + s2s3);
    double dotsum = 1 + s1s2 + s1s3 + s2s3;

    xt::xtensorf<double, xt::xshape<4>> s1xs3 = c_prod(s1, s3);
    double crosssum = xt::sum(s1xs3 * s2)[0];

    double ang = atan2(crosssum / rho, dotsum / rho);

    return ang;
}
