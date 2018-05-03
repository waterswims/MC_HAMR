#include "../includes/stdrand.hpp"

stdrand::std_d_unirand::std_d_unirand(int seed)
{
    distribution = std::uniform_real_distribution<double>(0.0, 1.0);
    this->change_seed(seed);
}

stdrand::std_i_unirand::std_i_unirand(int seed)
{
    distribution = std::uniform_int_distribution<int>(0, 1);
    this->change_seed(seed);
}

stdrand::std_normrand::std_normrand(double m, double sdin, int seed)
{
    distribution = std::normal_distribution<double>(m, sdin);
    this->change_seed(seed);
}

stdrand::std_lognormrand::std_lognormrand(double m, double sdin, int seed)
{
    distribution = std::lognormal_distribution<double>(m, sdin);
    this->change_seed(seed);
}
