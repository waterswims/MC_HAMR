#include "../includes/thermodynamics.hpp"
#include "../includes/functions.hpp"

#include <xtensor/xeval.hpp>
#include <xtensor/xview.hpp>

xt::xtensorf<double, xt::xshape<4>> Vsum, Vdiff, Vtemp, Vtemp2, Vloc_diff;

void particle::td::functionObject::setup(bool useJ, bool useD)
{
    std::function<void(field::field_type&, int, int)> insideE_func =
        [](field::field_type& lattice, int i, int j){};
    std::function<void(field::field_type&, int, int)> insidedE_func =
        [](field::field_type& lattice, int position, int j){};

    if(useJ)
    {
        insideE_func = [](field::field_type& lattice, int i, int j)
            {
                Vtemp += lattice.get_J(i, j) *
                    lattice.access(lattice.get_neigh(i)[j]);
            };
        insidedE_func = [](field::field_type& lattice,
            int position, int j)
            {
                Vsum += lattice.get_J(position, j) *
                        lattice.access(lattice.get_neigh(position)[j]);
            };
    }

    if(useD)
    {
        insideE_func = [insideE_func](field::field_type& lattice, int i, int j)
            {
                insideE_func(lattice, i, j);
                c_prod(lattice.get_D_vec(i, j),
                    lattice.access(lattice.get_neigh(i)[j]), Vtemp2);

                Vtemp += Vtemp2;
            };
        insidedE_func = [insidedE_func](field::field_type& lattice,
            int position, int j)
            {
                insidedE_func(lattice, position, j);
                c_prod(lattice.access(lattice.get_neigh(position)[j]),
                    lattice.get_D_vec(position, j), Vtemp);
                Vsum += Vtemp;
            };
    }

    dE_func = [insidedE_func](field::field_type& lattice, int position,
        xt::xtensorf<double, xt::xshape<4>>& H)
        {
            lattice.gen_rand();
            Vdiff = lattice.access(position) - lattice.get_rand();
            Vsum = {0, 0, 0, 0};

            for(int j = 0; j < lattice.get_neigh(position).size(); j++)
            {
                insidedE_func(lattice, position, j);
            }

            double dE = xt::eval(xt::sum((H + Vsum) * Vdiff))[0];

            return dE;
        };

    E_func = [insideE_func, this](field::field_type& lattice,
        xt::xtensorf<double, xt::xshape<4>>& H)
        {
            int Nspins = lattice.get_size();

            Vsum = {0, 0, 0, 0};

            for(int i = 0; i < Nspins; i++)
            {
                Vtemp = {0, 0, 0, 0};
                for(int j = 0; j < lattice.get_neigh(i).size(); j++)
                {
                    insideE_func(lattice, i, j);
                }

                Vsum += Vtemp * lattice.access(i);
            }

            Vtemp = calc_M(lattice);

            double E = xt::eval(xt::sum(-H * Vtemp - 0.5 * Vsum))[0];

            return E;
        };
}

xt::xtensorf<double, xt::xshape<4>> particle::td::functionObject::calc_M(
    particle::field::field_type& lattice)
{
    int Nspins = lattice.get_size();
    Vtemp = {0, 0, 0, 0};

    for(int i = 0; i < Nspins; i++)
    {
        Vtemp += lattice.access(i);
    }

    return Vtemp;
}

xt::xtensorf<double, xt::xshape<4>> particle::td::functionObject::calc_subM(
    particle::field::field_type& lattice,
    int subnumber)
{
    int Nspins = lattice.get_size();
    Vtemp = {0, 0, 0, 0};

    for(int i = 0; i < Nspins; i++)
    {
        int possum = xt::eval(xt::sum(lattice.get_loc(i)))[0];
        if (possum%2 == subnumber)
        {
            Vtemp += lattice.access(i);
        }
    }

    return Vtemp;
}

std::vector<double> particle::td::functionObject::calc_TC(
    particle::field::field_type& lattice)
{
    int Nspins = lattice.get_size();
    std::vector<double> out;
    int left, right, forward, backward;

    for(int i = 0; i < Nspins; i++)
    {
        left = -1;
        right = -1;
        forward = -1;
        backward = -1;

        if(out.size() < lattice.get_loc(i)[2]+1)
        {
            out.resize(lattice.get_loc(i)[2]+1);
        }
        for(int j = 0; j < Nspins; j++)
        {
            Vloc_diff = lattice.get_loc(j) - lattice.get_loc(i);
            if(Vloc_diff[2] == 0)
            {
                if (Vloc_diff[0] == 1 && Vloc_diff[1] == 0)
                {
                    left = j;
                }
                else if (Vloc_diff[0] == -1 && Vloc_diff[1] == 0)
                {
                    right = j;
                }
                else if (Vloc_diff[0] == 0 && Vloc_diff[1] == 1)
                {
                    forward = j;
                }
                else if (Vloc_diff[0] == 0 && Vloc_diff[1] == -1)
                {
                    backward = j;
                }
            }
        }

        if (right >= 0 && forward >= 0)
        {
            out[lattice.get_loc(i)[2]] += solid_angle(lattice.access(i),
                lattice.access(right), lattice.access(forward));
        }
        if (left >= 0 && backward >= 0)
        {
            out[lattice.get_loc(i)[2]] += solid_angle(lattice.access(i),
                lattice.access(left), lattice.access(backward));
        }
    }
    for(int i = 0; i < out.size(); i++)
    {
        out[i] /= 2 * M_PI;
    }

    return out;
}

double particle::td::solid_angle(const xt::xtensorf<double, xt::xshape<4>>& s1,
                const xt::xtensorf<double, xt::xshape<4>>& s2,
                const xt::xtensorf<double, xt::xshape<4>>& s3)
{
    double s1s2 = xt::eval(xt::sum(s1 * s2))[0];
    double s1s3 = xt::eval(xt::sum(s1 * s3))[0];
    double s2s3 = xt::eval(xt::sum(s2 * s3))[0];

    double rho = 2 * (1 + s1s2) * (1 + s1s3) * (1 + s2s3);
    double dotsum = 1 + s1s2 + s1s3 + s2s3;

    c_prod(s1, s3, Vtemp);
    double crosssum = xt::eval(xt::sum(Vtemp * s2))[0];

    double ang = atan2(crosssum / rho, dotsum / rho);

    return ang;
}
