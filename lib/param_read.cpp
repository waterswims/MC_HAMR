#include "../includes/param_read.hpp"
#include "../includes/array_alloc.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>

template <class T> T read_var(std::string v_name, std::string f_name)
{
    std::fstream file;
    file.open(f_name.c_str(), std::fstream::in);
    if(!file.is_open())
    {
        std::cout << "Input file not opened" << std::endl;
        exit(105);
    }
    std::string line;
    std::string test_v_name;
    T test_out;
    while(getline(file, line))
    {
        std::istringstream stream(line);
        stream >> test_v_name >> test_out;
        if (test_v_name == v_name)
        {
            file.close();
            return test_out;
        }
    }
    file.close();
    std::cout << "Variable name '" << v_name << "' not found" << std::endl;
    exit(105);
}

void read_all_vars(std::string f_name, double& size, double& J, double& k,
    bool& periodic, char& shape, char& hamil, int& Samp_steps, int& N_samp,
    int& Eq_steps, int& N_latts, double& beta, bool& distrib, double& amean,
    double& asd, std::string& temp_name, std::string& field_name, int& protocol,
    double& K, bool& print_latt)
{
    shape = read_var<char>("LATTSHAPE", f_name);
    switch(shape)
    {
        case 'w':
        case 'W':
        case 'x':
        case 'X':
            beta = read_var<double>("WEIBULLFACT", f_name);
            break;
        default:
            beta = 50;
    }
    hamil = read_var<char>("HAMILTONIAN", f_name);
    switch(hamil)
    {
        case 'i':
        case 'I':
        case 'h':
        case 'H':
            J = read_var<double>("EXCHANGE", f_name);
            break;
        case 's':
        case 'S':
            J = read_var<double>("EXCHANGE", f_name);
            K = read_var<double>("DMISTREN", f_name);
            break;
        default:
            J = 1;
            K = 0;
    }
    distrib = read_var<int>("ISDISTRIB", f_name);
    if(distrib)
    {
        amean = read_var<double>("MEANSIZE", f_name);
        asd = read_var<double>("SIZEDEV", f_name);
        size = amean;
    }
    else
    {
        size = read_var<double>("SIZE", f_name);
        amean = 0;
        asd = 0;
    }
    N_latts = read_var<int>("NLATTS", f_name);
    k = read_var<double>("BOLTZMANN", f_name);
    periodic = read_var<int>("ISPERIO", f_name);
    Samp_steps = read_var<int>("SAMPSWEEPS", f_name);
    Eq_steps = read_var<int>("EQSWEEPS", f_name);
    N_samp = read_var<int>("NSAMPS", f_name);
    temp_name = read_var<std::string>("TEMPNAME", f_name);
    field_name = read_var<std::string>("FIELDNAME", f_name);
    protocol = read_var<int>("PROTOCOL", f_name);
    print_latt = read_var<int>("PRINT_LATT", f_name);
}

void load_Hs_Ts(std::string Tname,
                float* &Ts,
                int& Tnum,
                std::string Hname,
                float* &Hs,
                int& Hnum)
{
    // load temps
    std::stringstream loadstream;
    std::string loadname;
    loadstream << "Temps/" << Tname << std::endl;
    loadstream >> loadname;

    std::ifstream f;
    f.open(loadname.c_str());
    double curr;
    bool cont = false;
    if(f >> curr) {cont = true;}
    int i = 0;
    for (; cont; i++)
    {
        if(f >> curr) {}
        else {cont = false;}
    }
    Tnum = i;
    f.close();

    Ts = alloc_1darr<float>(Tnum);

    f.open(loadname.c_str());
    cont = false;
    if(f >> curr) {cont = true;}
    for (int i = 0; cont; i++)
    {
        Ts[i] = curr;
        if(f >> curr) {}
        else {cont = false;}
    }
    f.close();

    // load fields
    loadstream << "Fields/" << Hname << std::endl;
    loadstream >> loadname;

    f.open(loadname.c_str());
    cont = false;
    if(f >> curr) {cont = true;}
    i = 0;
    for (; cont; i++)
    {
        if(f >> curr) {}
        else {cont = false;}
    }
    Hnum = i;
    f.close();

    Hs = alloc_1darr<float>(Hnum);

    f.open(loadname.c_str());
    cont = false;
    if(f >> curr) {cont = true;}
    for (int i = 0; cont; i++)
    {
        Hs[i] = curr;
        if(f >> curr) {}
        else {cont = false;}
    }
    f.close();
}
