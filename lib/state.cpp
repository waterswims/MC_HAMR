#include "../includes/state.hpp"
#include "../includes/td_funcs.hpp"

#ifdef __INTEL_COMPILER
#include "../includes/mklrand.hpp"
#define IRANDTYPE mklrand::mkl_irand
#define DRANDTYPE mklrand::mkl_drand
#else
#include "../includes/stdrand.hpp"
#define IRANDTYPE stdrand::std_i_unirand
#define DRANDTYPE stdrand::std_d_unirand
#endif

#include "../includes/functions.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <valarray>

extern IRANDTYPE st_rand_int;
extern DRANDTYPE st_rand_double;

state::state(double size, bool isPerio, bool ising, char shape_code, double J,
    double Hin, double k, double Temp, double K, double* args)
{
    s_code = shape_code;

    int H_size = 4;
    h_ind = 2;
    if(ising)
    {
        H_size = 1;
        h_ind = 0;
    }
    H.resize(H_size);
    H = 0;
    H[h_ind] = Hin;

    this->init_points(size, isPerio, J, K, args);

    k_b = k;
    this->change_temp(Temp);
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

void state::init_points(double size, bool isPerio, double J, double K, double* args)
{
    int pad = 1;
    double sizeab = size;
    double sizec = size;
    field = particle::field::field_type(false, isPerio, 3, size, J, K, "");
    switch (s_code)
    {
        case 's':
        case 'S':
            shape = new particle::shape::square;
            break;
        case 'w':
        case 'W':
            shape = new particle::shape::weibull((size), args[0]);
            break;
        case 'c':
        case 'C':
            shape = new particle::shape::cube;
            break;
        case 'x':
        case 'X':
            shape = new particle::shape::weibull((size), args[0]);
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
    std::valarray<int> posva(dim);
    posva = 0;
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
        dE = particle::funcs::calc_dE(field, choice, H);
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
    std::valarray<double> M = particle::funcs::calc_M(field);
    std::vector<double> M_out;
    for(int i=0; i < M.size(); i++)
    {
        M_out.push_back(M[i]);
    }
    return M_out;
}

std::vector<double> state::submag(int subnumber)
{
    std::valarray<double> M = particle::funcs::calc_subM(field, 0);
    std::vector<double> M_out;
    for(int i=0; i < M.size(); i++)
    {
        M_out.push_back(M[i]);
    }
    return M_out;
}

double state::energy()
{
    return particle::funcs::calc_E(field, H);
}

std::vector<double> state::tcharge()
{
    return particle::funcs::calc_TC(field);
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
    // std::cout << "T = " << T << ", k = " << k_b << ", beta = " << beta << std::endl;
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
