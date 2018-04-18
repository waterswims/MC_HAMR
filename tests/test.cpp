#include <gtest/gtest.h>

#include "mklrand_test.hpp"
#include "ising_test.hpp"
// #include "heis_test.hpp"
// #include "fept_test.hpp"
// #include "skyrm_test.hpp"

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
