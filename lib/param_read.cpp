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

void read_all_vars(std::string f_name, stateOptions& stOpt, simOptions& simOpt)
{
    stOpt.isIsing = read_var<int>("ISISING", f_name);

    stOpt.J = read_var<double>("EXCHANGE", f_name);
    stOpt.K = read_var<double>("DMISTREN", f_name);

    simOpt.distrib = read_var<int>("ISDISTRIB", f_name);
    if(simOpt.distrib)
    {
        simOpt.amean = read_var<double>("MEANSIZE", f_name);
        simOpt.asd = read_var<double>("SIZEDEV", f_name);
    }
    else
    {
        simOpt.amean = read_var<double>("SIZE", f_name);
        simOpt.asd = 0;
    }

    stOpt.shape_code = read_var<char>("LATTSHAPE", f_name);
    switch(stOpt.shape_code)
    {
        case 'w':
        case 'W':
        case 'x':
        case 'X':
            stOpt.beta = read_var<double>("WEIBULLFACT", f_name);
            stOpt.edgeSize = (simOpt.amean + 2*simOpt.asd) * 2 + 10;
            break;
        case 's':
        case 'S':
        case 'c':
        case 'C':
            stOpt.edgeSize = simOpt.amean + 2*simOpt.asd;
        break;
    }

    simOpt.N_latts = read_var<int>("NLATTS", f_name);
    stOpt.k = read_var<double>("BOLTZMANN", f_name);
    stOpt.isPerio = read_var<int>("ISPERIO", f_name);
    simOpt.Samp_steps = read_var<int>("SAMPSWEEPS", f_name);
    simOpt.Eq_steps = read_var<int>("EQSWEEPS", f_name);
    simOpt.N_samp = read_var<int>("NSAMPS", f_name);
    simOpt.tempFile = read_var<std::string>("TEMPNAME", f_name);
    simOpt.fieldFile = read_var<std::string>("FIELDNAME", f_name);
    simOpt.protocol = read_var<int>("PROTOCOL", f_name);
    simOpt.printLatt = read_var<int>("PRINT_LATT", f_name);
    stOpt.intFile = read_var<std::string>("INTERACTIONS", f_name);
}

void load_Hs_Ts(simOptions& simOpt,
                float* &Ts,
                int& Tnum,
                float* &Hs,
                int& Hnum)
{
    // load temps
    std::stringstream loadstream;
    std::string loadname;
    loadstream << "Temps/" << simOpt.tempFile << std::endl;
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
    loadstream << "Fields/" << simOpt.fieldFile << std::endl;
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
