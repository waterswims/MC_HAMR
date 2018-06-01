#ifndef _TESTFUNCS
#define _TESTFUNCS

#include "../includes/field_type.hpp"
#include "../includes/stdrand.hpp"
#define IRANDTYPE stdrand::std_i_unirand
#define DRANDTYPE stdrand::std_d_unirand
#define LNRANDTYPE stdrand::std_lognormrand

#include <vector>
#include <iostream>
#include <cmath>
#include <xtensor/xeval.hpp>

IRANDTYPE st_rand_int(1e5, 1);
DRANDTYPE st_rand_double(1e5, 2);
LNRANDTYPE rand_ln(0, 0.25, 1e5, 3);

particle::field::field_type gen_fm(int d, bool ising, double J, double D, bool perio)
{
    particle::field::field_type
        field(ising, perio, d, 10, J, D, "Js/3d_nn_bloch.dat");
    xt::xtensorf<int, xt::xshape<4>> pos;
    pos = {0, 0, 0, 0};
    for(int i = 0; i < pow(10, d); i++)
    {
        field.add_spin(pos);
        pos[d-1]++;
        for(int j=d-1; j > 0; j--)
        {
            if(pos[j] == 10)
            {
                pos[j] = 0;
                pos[j-1]++;
            }
        }
    }

    field.set_neigh();

    return field;
}

particle::field::field_type gen_afm(int d, bool ising, double J, double D, bool perio)
{
    particle::field::field_type
        field(ising, false, d, 10, J, D, "Js/3d_nn_bloch.dat");
    xt::xtensorf<int, xt::xshape<4>> pos;
    pos = {0, 0, 0, 0};
    for(int i = 0; i < pow(10, d); i++)
    {
        field.add_spin(pos);
        if(xt::eval(xt::sum(pos))[0]%2 == 1)
        {
            field.set_down(i);
        }
        pos[d-1]++;
        for(int j=d-1; j > 0; j--)
        {
            if(pos[j] == 10)
            {
                pos[j] = 0;
                pos[j-1]++;
            }
        }
    }

    field.set_neigh();

    return field;
}

particle::field::field_type gen_skyrm(int d, double J, double D)
{
    int edgesize = 30;
    particle::field::field_type
        field(false, false, d, edgesize, J, D, "Js/3d_nn_bloch.dat");
    xt::xtensorf<int, xt::xshape<4>> pos;
    pos = {0, 0, 0, 0};
    xt::xtensorf<double, xt::xshape<4>> spin;
    spin = {0, 0, 0, 0};
    int mid = edgesize / 2;
    float sk_R = mid - 3;
    float k = M_PI / sk_R;
    for(int i = 0; i < pow(edgesize, d); i++)
    {
        field.add_spin(pos);

        float dx = pos[0] - mid;
        float dy = pos[1] - mid;
        float r2 = pow(dx, 2) + pow(dy, 2);
        float r = pow(r2, 0.5);
        if(r >= sk_R)
        {
            spin = {0, 0, 0, 0};
            field.set_spin(i, spin);
        }
        else
        {
            float phi = atan2(dx, dy);
            spin[0] = -sin(k*r)*cos(phi);
            spin[1] = sin(k*r)*sin(phi);
            spin[2] = cos(k*r);
            field.set_spin(i, spin);

        }

        pos[d-1]++;
        for(int j=d-1; j > 0; j--)
        {
            if(pos[j] == edgesize)
            {
                pos[j] = 0;
                pos[j-1]++;
            }
        }
    }

    field.set_neigh();

    return field;
}

#endif
