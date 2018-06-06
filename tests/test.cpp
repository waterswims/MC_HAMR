#include <gtest/gtest.h>

#include "test_functions.hpp"
#include "rand_test.hpp"
#include "ising_test.hpp"
#include "heis_test.hpp"
// #include "fept_test.hpp"

int main(int argc, char **argv)
{
    std::cout << "Generating fields" << std::endl;
    isingFMField = gen_fm(3, true, 1, 0, false);
    isingFMFieldPerio = gen_fm(3, true, 1, 0, true);
    isingAFMField = gen_afm(3, true, 1, 0, false);
    heisFMField = gen_fm(3, false, 1, 0, false);
    heisAFMField = gen_afm(3, false, 1, 0, false);
    skyrmField = gen_skyrm(3, 1, 0.75);
    exchangeOnly.setup(true, false);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
