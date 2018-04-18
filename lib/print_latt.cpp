#include "../includes/print_latt.hpp"
#include <cstdlib>
#include <iostream>
#include <sstream>

std::string av_latt_name(const int protocol,
                    const int var1,
                    const int var2)
{
    int H, T;
    switch(protocol)
    {
        case 1:
        H = var1;
        T = var2;
        break;
        case 2:
        H = var2;
        T = var1;
        break;
    }

    std::stringstream nstream;
    std::string name;

    nstream << "/Av_Latt/T_" << T << "-H_" << H;
    nstream >> name;

    return name;
}

std::string sing_latt_name(const int protocol,
                    const int var1,
                    const int var2)
{
    int H, T;
    switch(protocol)
    {
        case 1:
        H = var1;
        T = var2;
        break;
        case 2:
        case 4:
        H = var2;
        T = var1;
        break;
    }

    std::stringstream nstream;
    std::string name;

    nstream << "/Sing_Latt/T_" << T << "-H_" << H;
    nstream >> name;

    return name;
}
