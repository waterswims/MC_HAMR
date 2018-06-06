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
            Vdiff = lattice.access(position) - lattice.get_rand();
            Vsum = {0, 0, 0, 0};

            for(int j = 0; j < lattice.get_neigh(position).size(); j++)
            {
                insidedE_func(lattice, position, j);
            }

            auto res = (H + Vsum) * Vdiff;

            double dE = res[0] + res[1] + res[2];

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

            auto res = -H * Vtemp - 0.5 * Vsum;

            double E = res[0] + res[1] + res[2];

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
    int d = lattice.get_dim();
    int tcsize = 1;
    if(d == 3)
    {
        tcsize = lattice.get_edge();
    }
    std::vector<double> out(tcsize);

    for(int i = 0; i < Nspins; i++)
    {
        if (lattice.get_adj(i)[1] >= 0 && lattice.get_adj(i)[2] >= 0)
        {
            out[lattice.get_loc(i)[2]] += solid_angle(lattice.access(i),
                lattice.access(lattice.get_adj(i)[1]),
                lattice.access(lattice.get_adj(i)[2]));
        }
        if (lattice.get_adj(i)[0] >= 0 && lattice.get_adj(i)[3] >= 0)
        {
            out[lattice.get_loc(i)[2]] += solid_angle(lattice.access(i),
                lattice.access(lattice.get_adj(i)[0]),
                lattice.access(lattice.get_adj(i)[3]));
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
    auto tempres = s1 * s2;
    double s1s2 = tempres[0] + tempres[1] + tempres[2];
    auto tempres2 = s1 * s3;
    double s1s3 = tempres2[0] + tempres2[1] + tempres2[2];
    auto tempres3 = s2 * s3;
    double s2s3 = tempres3[0] + tempres3[1] + tempres3[2];

    double rho = 2 * (1 + s1s2) * (1 + s1s3) * (1 + s2s3);
    double dotsum = 1 + s1s2 + s1s3 + s2s3;

    c_prod(s1, s3, Vtemp);
    auto tempres4 = Vtemp * s2;
    double crosssum = tempres4[0] + tempres4[1] + tempres4[2];

    double ang = atan2(crosssum / rho, dotsum / rho);

    return ang;
}
