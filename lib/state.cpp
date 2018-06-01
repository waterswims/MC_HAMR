#include "../includes/state.hpp"
#include "../includes/stdrand.hpp"
#define IRANDTYPE stdrand::std_i_unirand
#define DRANDTYPE stdrand::std_d_unirand

#include "../includes/functions.hpp"
#include <iostream>
#include <fstream>
#include <cmath>

extern IRANDTYPE st_rand_int;
extern DRANDTYPE st_rand_double;

state::state(stateOptions opt)
{
    td_funcs.setup((opt.J!=0), (opt.K!=0));

    s_code = opt.shape_code;

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
    beta = other.beta;
    k_b = other.k_b;
    num = other.num;
    s_code = other.s_code;
    this->copy_points(other);
}

state& state::operator=(const state& other)
{
    beta = other.beta;
    k_b = other.k_b;
    num = other.num;
    s_code = other.s_code;
    this->copy_points(other);
}

void state::init_points(stateOptions opt)
{
    int pad = 1;
    double sizeab = opt.size;
    double sizec = opt.size;
    int size;
    switch (s_code)
    {
        case 's':
        case 'S':
            shape = new particle::shape::square;
            field = particle::field::field_type(opt.isIsing,
                opt.isPerio, 2, opt.size, opt.J, opt.K, opt.intFile);
            break;
        case 'w':
        case 'W':
            size = opt.size * 2 + 10;
            shape = new particle::shape::weibull((size), opt.beta);
            field = particle::field::field_type(opt.isIsing,
                opt.isPerio, 2, size, opt.J, opt.K, opt.intFile);
            break;
        case 'c':
        case 'C':
            shape = new particle::shape::cube;
            field = particle::field::field_type(opt.isIsing,
                opt.isPerio, 3, opt.size, opt.J, opt.K, opt.intFile);
            break;
        case 'x':
        case 'X':
            size = opt.size * 2 + 10;
            shape = new particle::shape::weibull((size), opt.beta);
            field = particle::field::field_type(opt.isIsing,
                opt.isPerio, 3, size, opt.J, opt.K, opt.intFile);
            break;
        default:
            std::cerr << "Incorrect shape code, exiting" << std::endl;
            exit(103);
    }
}

void state::copy_points(const state& other)
{
    field = other.field;
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
    int size = field.get_size();
    int dim = field.get_dim();
    std::vector<int> pos(dim, 0);
    xt::xtensorf<int, xt::xshape<4>> posva = {0, 0, 0, 0};
    for(int i = 0; i < pow(size, dim); i++)
    {
        bool fillspin = shape->check(pos, size);
        int possum = sum(pos);
        if (fillspin)
        {
            field.add_spin(posva);
            num++;
            if (possum%2 == 0){snum++;}
        }
        for(int j; j < dim-1; j++)
        {
            pos[dim-j]++;
            pos[dim-j]++;
            if(pos[dim-j] == 10)
            {
                pos[dim-j] = 0;
                pos[dim-j-1]++;
                posva[dim-j] = 0;
                posva[dim-j-1]++;
            }
        }
    }

    field.set_neigh();
}

void state::equil(int iter)
{
    int dim = field.get_dim();
    int size = field.get_size();
    int choice;

    // create some variables
    double dE = 0;
    double log_eta = 0;

    for (int i=0; i<iter; i++)
    {
        choice = int(st_rand_double.gen() * size);

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
    int dim = field.get_dim();

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
