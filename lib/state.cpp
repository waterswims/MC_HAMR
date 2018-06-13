#include "../includes/state.hpp"
#include "../includes/stdrand.hpp"
#define IRANDTYPE stdrand::std_i_unirand
#define DRANDTYPE stdrand::std_d_unirand

#include "../includes/functions.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <xtensor/xio.hpp>

extern IRANDTYPE st_rand_int;
extern DRANDTYPE st_rand_double;

state::state(stateOptions opt)
{
    td_funcs.setup((opt.J!=0), (opt.K!=0));

    s_code = opt.shape_code;
    edgesize = opt.edgeSize;

    h_ind = 2;
    if(opt.isIsing)
    {
        h_ind = 0;
    }
    H = {0, 0, 0, 0};
    H[h_ind] = opt.H;

    this->init_points(opt);

    k_b = opt.k;
    this->change_temp(opt.T);
    this->init_lattice();
}

state::state(const state& other)
{
    this->copy_points(other);
}

state& state::operator=(const state& other)
{
    this->copy_points(other);
    return *this;
}

void state::init_points(stateOptions opt)
{
    int pad = 1;
    double sizeab = opt.size;
    double sizec = opt.size;
    int d;
    switch (s_code)
    {
        case 's':
        case 'S':
            d = 2;
            shape = new particle::shape::square;
            break;
        case 'w':
        case 'W':
            d = 2;
            shape = new particle::shape::weibull(opt.size, opt.beta);
            break;
        case 'c':
        case 'C':
            d = 3;
            shape = new particle::shape::cube;
            break;
        case 'x':
        case 'X':
            d = 3;
            shape = new particle::shape::weibull(opt.size, opt.beta);
            break;
        default:
            std::cerr << "Incorrect shape code, exiting" << std::endl;
            exit(103);
    }
    field = particle::field::field_type(opt.isIsing, opt.isPerio, d,
        opt.edgeSize, opt.J, opt.K, opt.intFile);
}

void state::copy_points(const state& other)
{
    field = other.field;
    td_funcs = other.td_funcs;
    beta = other.beta;
    k_b = other.k_b;
    num = other.num;
    s_code = other.s_code;
    edgesize = other.edgesize;
    H = other.H;
    h_ind = other.h_ind;
    switch (s_code)
    {
        case 's':
        case 'S':
            shape = new particle::shape::square;
            break;
        case 'w':
        case 'W':
            shape = new particle::shape::weibull(*(other.shape));
            break;
        case 'c':
        case 'C':
            shape = new particle::shape::cube;
            break;
        case 'x':
        case 'X':
            shape = new particle::shape::weibull(*(other.shape));
            break;
        default:
            std::cerr << "Incorrect shape code, exiting" << std::endl;
            exit(103);
    }
}

state::~state()
{
    delete shape;
}

void state::init_lattice()
{
    num = 0;
    snum = 0;
    s4num = 0;
    int dim = field.get_dim();
    std::vector<int> pos(dim, 0);
    xt::xtensorf<int, xt::xshape<4>> posva = {0, 0, 0, 0};
    for(int i = 0; i < pow(edgesize, dim); i++)
    {
        bool fillspin = shape->check(pos, edgesize);
        int possum = sum(pos);
        if (fillspin)
        {
            field.add_spin(posva);
            num++;

            if (possum%2 == 0){snum++;}

            int possum2 = posva[0] + posva[1];
            int posdiff = -posva[0] + posva[1] + (edgesize - edgesize%4);
            if (possum2%4 == 0 && posdiff%4 == 0 && posva[2]%4 == 0)
            {
                s4num++;
            }
            else if (possum2%4 == 2 && posdiff%4 == 2 && posva[2]%4 == 2)
            {
                s4num++;
            }
        }
        pos[dim-1]++;
        posva[dim-1]++;
        for(int j=dim-1; j > 0; j--)
        {
            if(pos[dim-j] == edgesize)
            {
                pos[dim-j] = 0;
                pos[dim-j-1]++;
                posva[dim-j] = 0;
                posva[dim-j-1]++;
            }
        }
    }
    field.all_rand();
    field.set_neigh();
}

void state::equil(int iter)
{
    int size = field.get_size();
    int choice;

    // create some variables
    double dE = 0, spare;
    double log_eta = 0;

    for (int i=0; i<iter; i++)
    {
        choice = int(st_rand_double.gen() * size);

        field.gen_rand();

        //check dE
        dE = td_funcs.calc_dE(field, choice, H);

        //check if flip
        if(dE <= 0)
        {
            field.set_rand(choice);
        }
        else
        {
            log_eta = log(st_rand_double.gen());
			if ((-dE * beta) > log_eta)
			{
                field.set_rand(choice);
			}
        }
    }
}

std::vector<double> state::magnetisation()
{
    xt::xtensorf<double, xt::xshape<4>> M = td_funcs.calc_M(field);
    std::vector<double> M_out;
    for(int i=0; i < 4; i++)
    {
        M_out.push_back(M[i]);
    }
    return M_out;
}

std::vector<double> state::submag(int subnumber)
{
    xt::xtensorf<double, xt::xshape<4>> M = td_funcs.calc_subM(field, 0);
    std::vector<double> M_out;
    for(int i=0; i < 4; i++)
    {
        M_out.push_back(M[i]);
    }
    return M_out;
}

std::vector<double> state::sub4mag()
{
    xt::xtensorf<double, xt::xshape<4>> M = td_funcs.calc_sub4M(field);
    std::vector<double> M_out;
    for(int i=0; i < 4; i++)
    {
        M_out.push_back(M[i]);
    }
    return M_out;
}

double state::energy()
{
    return td_funcs.calc_E(field, H);
}

std::vector<double> state::tcharge()
{
    return td_funcs.calc_TC(field);
}

int state::num_spins()
{
    return num;
}

int state::sub_num(int subnumber)
{
    return snum;
}

int state::sub4_num()
{
    return s4num;
}

void state::change_temp(double T)
{
    if (T <= 0)
    {
        std::cerr << "Invalid temperature, exiting" << std::endl;
        exit(104);
    }
    beta = 1.0 / (k_b * T);
}

void state::change_field(double Hin)
{
    H[h_ind] = Hin;
}

void state::print_latt()
{

}

void state::ptf(std::string fname, std::string arrname)
{
    field.print(fname, arrname);
}

void state::add_to_av(particle::field::field_type& other_field)
{
    int size = field.get_size();

    for(int i = 0; i < size; i++)
    {
        other_field.access(i) += field.access(i);
    }
}

void state::change_v1(int protocol, double v1)
{
    switch(protocol)
    {
        case 1:
        case 3:
        this->change_field(v1);
        break;
        case 2:
        case 4:
        this->change_temp(v1);
        break;
    }
}

void state::change_v2(int protocol, double v2)
{
    switch(protocol)
    {
        case 1:
        case 3:
        this->change_temp(v2);
        break;
        case 2:
        case 4:
        this->change_field(v2);
        break;
    }
}

void state::send_latt_data(int dest_rank)
{
    field.send_data(dest_rank);
}

void state::recv_latt_data(int src_rank)
{
    field.recv_data(src_rank);
}
