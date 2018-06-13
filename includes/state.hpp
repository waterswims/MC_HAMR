#ifndef _STATE
#define _STATE

#include "param_read.hpp"
#include "field_type.hpp"
#include "thermodynamics.hpp"
#include "shape.hpp"

#include <string.h>
#include <vector>

class state
{
private:
    double beta, k_b;
    xt::xtensorf<double, xt::xshape<4>> H;
    int num, snum, s4num, h_ind, edgesize;
    char s_code;
    particle::field::field_type field;
    particle::td::functionObject td_funcs;
    particle::shape::shape_type* shape;

public:
    state(){}
    state(stateOptions opt);
    state(const state& other);
    ~state();
    particle::field::field_type get_field() {return field;}
    void copy_points(const state& other);
    void init_points(stateOptions opt);
    void equil(int iter);
    std::vector<double> magnetisation();
    std::vector<double> submag(int subnumber);
    std::vector<double> sub4mag();
    double energy();
    std::vector<double> tcharge();
    int num_spins();
    int get_size() {return edgesize;}
    int sub_num(int subnumber);
    int sub4_num();
    void init_lattice();
    void change_temp(double T);
    void change_field(double Hin);
    state& operator=(const state& other);
    void print_latt();
    void ptf(std::string fname, std::string arrname);
    void add_to_av(particle::field::field_type& other_field);
    void change_v1(int protocol, double v1);
    void change_v2(int protocol, double v2);
    void send_latt_data(int dest_rank);
    void recv_latt_data(int src_rank);
};

#endif
