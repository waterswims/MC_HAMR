#include "../includes/stdrand.hpp"

stdrand::std_d_unirand::std_d_unirand(int seed) : std_randbase()
{
    this->change_seed(seed);
}

stdrand::std_i_unirand::std_i_unirand(int seed) : std_randbase()
{
    this->change_seed(seed);
}

stdrand::std_normrand::std_normrand(double m, double sdin, int seed) : std_randbase()
{
    M = m;
    SD = sdin;
    this->change_seed(seed);
}

stdrand::std_lognormrand::std_lognormrand(double m, double sdin, int seed) : std_randbase()
{
    M = m;
    SD = sdin;
    this->change_seed(seed);
}

stdrand::std_sphere::std_sphere(int seed) : std_randbase()
{
    this->change_seed(seed);
}
