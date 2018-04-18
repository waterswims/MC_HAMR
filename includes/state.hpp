#ifndef _STATE
#define _STATE

#include "field_type.hpp"
#include "shape.hpp"

#include <string.h>
#include <vector>
#include <valarray>

class state
{
private:
    double beta, k_b;
    std::valarray<double> H;
    int num, snum, h_ind;
    char s_code;
    particle::field::field_type field;
    particle::shape::shape_type* shape;

public:
    state(){}
    state(double size, bool isPerio, bool ising, char shape_code, double J,
        double Hin, double k, double Temp, double K, double* args);
    state(const state& other);
    ~state();
    void copy_points(const state& other);
    void init_points(double size, bool isPerio, double J, double K, double* args);
    void equil(int iter);
    std::vector<double> magnetisation();
    std::vector<double> submag(int subnumber);
    double energy();
    std::vector<double> tcharge();
    int num_spins();
    int sub_num(int subnumber);
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
