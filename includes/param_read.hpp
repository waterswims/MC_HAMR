#include <string>

template <class T> T read_var(std::string v_name, std::string f_name);

void read_all_vars(std::string f_name, double& size, double& J, double& k,
    bool& periodic, char& shape, char& hamil, int& Samp_steps, int& N_samp,
    int& Eq_steps, int& N_latts, double& beta, bool& distrib, double& amean,
    double& asd, std::string& temp_name, std::string& field_name, int& protocol,
    double& K, bool& print_latt);

void load_Hs_Ts(std::string Tname,
                float* &Ts,
                int& Tnum,
                std::string Hname,
                float* &Hs,
                int& Hnum);
