#ifndef _PARAM
#define _PARAM

#include <string>

struct stateOptions
{
    double size;
    double edgeSize;
    bool isPerio;
    bool isIsing;
    char shape_code;
    double J;
    double K;
    double H;
    double k;
    double T;
    double beta;
    std::string intFile;
};

struct simOptions
{
    std::string tempFile;
    std::string fieldFile;
    std::string outFile;
    int Samp_steps;
    int N_samp;
    int Eq_steps;
    int N_latts;
    int protocol;
    bool printLatt;
    bool distrib;
    double amean;
    double asd;
    double lmean;
    double lsd;
};

template <class T> T read_var(std::string v_name, std::string f_name);

void read_all_vars(std::string f_name, stateOptions& stOpt, simOptions& simOpt);

void load_Hs_Ts(simOptions& simOpt,
                float* &Ts,
                int& Tnum,
                float* &Hs,
                int& Hnum);

#endif
