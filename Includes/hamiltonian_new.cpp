#include "hamiltonian_new.hpp"
#include "array_alloc.hpp"
#include "functions.h"

#include <iostream>
#include <cstdlib>

///////////////////////////////////
// Ising model hamiltonian
///////////////////////////////////

ham_ising::ham_ising(ham_type& other)
{
    H = other.get_J();
    J = other.get_H();
}

ham_ising& ham_ising::operator=(ham_type& other)
{
    H = other.get_J();
    J = other.get_H();
    return *this;
}

double ham_ising::calc_E(field_type* lattice)
{
    int sum=0;
    int start = 0;
    if(lattice->get_perio())
    {
        start++;
    }
    int dim = lattice->get_dim();
    vector<int> pos(dim);
    for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        *it = start;
    }
    bool finished = false;
    int curr;
    int H_sum(0), J_sum(0);
    adj = alloc_1darr<int>(dim*2);
    while (!finished)
    {
        lattice->i_adjacent(pos, adj);
        lattice->i_next(finished, pos, curr);
        H_sum += curr;
        for (int i = 0; i < 2*dim; i++)
        {
            J_sum += curr*adj[i];
        }
    }
    dealloc_1darr<int>(dim*2, adj);

    double E = -H * H_sum - 0.5 * J * J_sum;
    return E;
}

vector<double> ham_ising::calc_M(field_type* lattice)
{
    int sum=0;
    int start = 0;
    if(lattice->get_perio())
    {
        start++;
    }
    int dim = lattice->get_dim();
    vector<int> pos(dim);
    for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        *it = start;
    }
    bool finished = false;
    int curr;
    while (!finished)
    {
        lattice->i_next(finished, pos, curr);
        sum += curr;
    }
    vector<double> mag(1);
    mag[0] = sum;
    return mag;
}

vector<double> ham_ising::calc_subM(field_type* lattice, int subnumber)
{
    if (subnumber > 1)
    {
        cerr << "subnumber is incorrect" << endl;
        exit(304);
    }
    int tsum=0;
    int start = 0;
    if(lattice->get_perio())
    {
        start++;
    }
    int dim = lattice->get_dim();
    vector<int> pos(dim);
    for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        *it = start;
    }
    bool finished = false;
    int possum = 0;
    int curr;
    while (!finished)
    {
        possum = sum(pos);
        lattice->i_next(finished, pos, curr);
        if (possum%2 == subnumber)
        {
            tsum += curr;
        }
    }

    vector<double> mag(1);
    mag[0] = tsum;
    return mag;
}

double ham_ising::dE(field_type* lattice, vector<int>& position)
{
    int val;
    lattice->i_access(position, val);
    double dEH = 2 * H * val;
    int sum = 0;
    int dim = lattice->get_dim();
    adj = alloc_1darr<int>(dim*2);
    lattice->i_adjacent(position, adj);
    for(int i = 0; i < dim*2; i++)
    {
        sum += adj[i];
    }
    dealloc_1darr<int>(dim*2, adj);

    double dE = dEH + 2 * J * val * sum;
    return dE;
}

///////////////////////////////////
// Heis model hamiltonian
///////////////////////////////////

ham_heis::ham_heis(double Hin, double Jin)
{
    H.resize(3);
    J.resize(3);
    H[0] = 0;
    H[1] = 0;
    H[2] = Hin;
    J[0] = Jin;
    J[1] = Jin;
    J[2] = Jin;
    vsum.resize(3);
    curr.resize(3);
    H_sum.resize(3);
    J_sum.resize(3);
    test.resize(3);
}

ham_heis::ham_heis(ham_type& other)
{
    H = other.get_Hs();
    J = other.get_Js();
    vsum.resize(3);
    curr.resize(3);
    H_sum.resize(3);
    J_sum.resize(3);
    test.resize(3);
}

ham_heis& ham_heis::operator=(ham_type& other)
{
    H = other.get_Hs();
    J = other.get_Js();
    vsum.resize(3);
    curr.resize(3);
    H_sum.resize(3);
    J_sum.resize(3);
    test.resize(3);

    return *this;
}

double ham_heis::calc_E(field_type* lattice)
{
    int sum=0;
    int start = 0;
    if(lattice->get_perio())
    {
        start++;
    }
    int dim = lattice->get_dim();
    pos.resize(dim);
    for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        *it = start;
    }
    bool finished = false;
    H_sum[0] = 0;
    H_sum[1] = 0;
    H_sum[2] = 0;
    J_sum[0] = 0;
    J_sum[1] = 0;
    J_sum[2] = 0;
    int arrsize = dim*2;
    adj = alloc_2darr<double>(3, arrsize);
    while (!finished)
    {
        lattice->h_adjacent(pos, adj);
        lattice->h_next(finished, pos, curr);
        for (int j = 0; j < 3; j++)
        {
            H_sum[j] += curr[j];
            for (int i = 0; i < arrsize; i++)
            {
                J_sum[j] += curr[j]*adj[j][i];
            }
        }
    }
    dealloc_2darr<double>(3, arrsize, adj);
    for (int j = 0; j < 3; j++)
    {
        H_sum[j] = H_sum[j] * H[j];
        J_sum[j] = J_sum[j] * J[j];
    }
    double E = -(H_sum[0] + H_sum[1] + H_sum[2]) - 0.5*(J_sum[0] + J_sum[1] + J_sum[2]);

    return E;
}

vector<double> ham_heis::calc_M(field_type* lattice)
{
    vsum[0] = 0;
    vsum[1] = 0;
    vsum[2] = 0;
    int start = 0;
    if(lattice->get_perio())
    {
        start++;
    }
    int dim = lattice->get_dim();
    pos.resize(dim);
    for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        *it = start;
    }
    bool finished = false;
    while (!finished)
    {
        lattice->h_next(finished, pos, curr);
        vsum[0] += curr[0];
        vsum[1] += curr[1];
        vsum[2] += curr[2];
    }
    return vsum;
}

vector<double> ham_heis::calc_subM(field_type* lattice, int subnumber)
{
    vsum[0] = 0;
    vsum[1] = 0;
    vsum[2] = 0;
    int start = 0;
    if(lattice->get_perio())
    {
        start++;
    }
    int dim = lattice->get_dim();
    pos.resize(dim);
    for (vector<int>::iterator it = pos.begin(); it != pos.end(); it++)
    {
        *it = start;
    }
    bool finished = false;
    int possum = 0;
    while (!finished)
    {
        possum = sum(pos);
        lattice->h_next(finished, pos, curr);
        if (possum%2 == subnumber)
        {
            vsum[0] += curr[0];
            vsum[1] += curr[1];
            vsum[2] += curr[2];
        }
    }
    return vsum;
}

double ham_heis::dE(field_type* lattice, vector<int>& position)
{
    rand_spin_h(test[0], test[1], test[2]);
    lattice->h_access(position, curr);
    double cmp1 = curr[0] - test[0];
    double cmp2 = curr[1] - test[2];
    double cmp3 = curr[2] - test[3];
    double dEH = H[0] * cmp1 + H[1] * cmp2 + H[2] * cmp3;

    int dim = lattice->get_dim();
    int arrsize = dim*2;
    adj = alloc_2darr<double>(3, arrsize);
    lattice->h_adjacent(position, adj);
    for(int j = 0; j < 3; j++)
    {
        vsum[j] = 0;
        for(int i = 0; i < arrsize; i++)
        {
            vsum[j] += adj[j][i];
        }
    }
    dealloc_2darr<double>(3, arrsize, adj);

    double dE = dEH + (J[0] * cmp1 * vsum[0] + J[1] * cmp2 * vsum[1] + J[2] * cmp3 * vsum[2]);

    return dE;
}
