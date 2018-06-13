#ifndef _HEISTEST
#define _HEISTEST

#include "test_functions.hpp"

xt::xtensorf<double, xt::xshape<4>> heisH = {0, 0, 0.1, 0};
particle::field::field_type heisFMField;
particle::field::field_type heisAFMField;
particle::field::field_type skyrmField;

///////////////////////////////////////////////////////
// Heisenberg model tests - 2D
///////////////////////////////////////////////////////

// TEST(Heis_model, 2d_energy_zero_field)
// {
//     field_2d_h field = gen_2d_heis_fm(1, 0, 0);
//     ham_heis hamil(0, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(-180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_fm(-1, 0, 0);
//     EXPECT_EQ(-180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_fm(0, 1, 0);
//     EXPECT_EQ(-180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_fm(0, -1, 0);
//     EXPECT_EQ(-180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_fm(0, 0, 1);
//     EXPECT_EQ(-180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_fm(0, 0, -1);
//     EXPECT_EQ(-180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_afm(1, 0, 0);
//     EXPECT_EQ(180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_afm(-1, 0, 0);
//     EXPECT_EQ(180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_afm(0, 1, 0);
//     EXPECT_EQ(180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_afm(0, -1, 0);
//     EXPECT_EQ(180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_afm(0, 0, 1);
//     EXPECT_EQ(180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_afm(0, 0, -1);
//     EXPECT_EQ(180, hamil.calc_E(&field));
// }
//
// TEST(Heis_model, 2d_energy_ext_field)
// {
//     field_2d_h field = gen_2d_heis_fm(1, 0, 0);
//     ham_heis hamil(0.1, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(-180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_fm(-1, 0, 0);
//     EXPECT_EQ(-180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_fm(0, 1, 0);
//     EXPECT_EQ(-180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_fm(0, -1, 0);
//     EXPECT_EQ(-180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_fm(0, 0, 1);
//     EXPECT_EQ(-190, hamil.calc_E(&field));
//
//     field = gen_2d_heis_fm(0, 0, -1);
//     EXPECT_EQ(-170, hamil.calc_E(&field));
//
//     field = gen_2d_heis_afm(1, 0, 0);
//     EXPECT_EQ(180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_afm(-1, 0, 0);
//     EXPECT_EQ(180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_afm(0, 1, 0);
//     EXPECT_EQ(180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_afm(0, -1, 0);
//     EXPECT_EQ(180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_afm(0, 0, 1);
//     EXPECT_EQ(180, hamil.calc_E(&field));
//
//     field = gen_2d_heis_afm(0, 0, -1);
//     EXPECT_EQ(180, hamil.calc_E(&field));
// }
//
// TEST(Heis_model, 2d_mag)
// {
//     field_2d_h field = gen_2d_heis_fm(1, 0, 0);
//     ham_heis hamil(0, 1);
//     hamil.init_dim(&field);
//     std::vector<double> mag = hamil.calc_M(&field);
//     EXPECT_EQ(100, mag[0]);
//     EXPECT_EQ(0, mag[1]);
//     EXPECT_EQ(0, mag[2]);
//
//     field = gen_2d_heis_fm(-1, 0, 0);
//     mag = hamil.calc_M(&field);
//     EXPECT_EQ(-100, mag[0]);
//     EXPECT_EQ(0, mag[1]);
//     EXPECT_EQ(0, mag[2]);
//
//     field = gen_2d_heis_fm(0, 1, 0);
//     mag = hamil.calc_M(&field);
//     EXPECT_EQ(0, mag[0]);
//     EXPECT_EQ(100, mag[1]);
//     EXPECT_EQ(0, mag[2]);
//
//     field = gen_2d_heis_fm(0, -1, 0);
//     mag = hamil.calc_M(&field);
//     EXPECT_EQ(0, mag[0]);
//     EXPECT_EQ(-100, mag[1]);
//     EXPECT_EQ(0, mag[2]);
//
//     field = gen_2d_heis_fm(0, 0, 1);
//     mag = hamil.calc_M(&field);
//     EXPECT_EQ(0, mag[0]);
//     EXPECT_EQ(0, mag[1]);
//     EXPECT_EQ(100, mag[2]);
//
//     field = gen_2d_heis_fm(0, 0, -1);
//     mag = hamil.calc_M(&field);
//     EXPECT_EQ(0, mag[0]);
//     EXPECT_EQ(0, mag[1]);
//     EXPECT_EQ(-100, mag[2]);
//
//     field = gen_2d_heis_afm(1, 0, 0);
//     mag = hamil.calc_M(&field);
//     EXPECT_EQ(0, mag[0]);
//     EXPECT_EQ(0, mag[1]);
//     EXPECT_EQ(0, mag[2]);
//
//     field = gen_2d_heis_afm(-1, 0, 0);
//     mag = hamil.calc_M(&field);
//     EXPECT_EQ(0, mag[0]);
//     EXPECT_EQ(0, mag[1]);
//     EXPECT_EQ(0, mag[2]);
//
//     field = gen_2d_heis_afm(0, 1, 0);
//     mag = hamil.calc_M(&field);
//     EXPECT_EQ(0, mag[0]);
//     EXPECT_EQ(0, mag[1]);
//     EXPECT_EQ(0, mag[2]);
//
//     field = gen_2d_heis_afm(0, -1, 0);
//     mag = hamil.calc_M(&field);
//     EXPECT_EQ(0, mag[0]);
//     EXPECT_EQ(0, mag[1]);
//     EXPECT_EQ(0, mag[2]);
//
//     field = gen_2d_heis_afm(0, 0, 1);
//     mag = hamil.calc_M(&field);
//     EXPECT_EQ(0, mag[0]);
//     EXPECT_EQ(0, mag[1]);
//     EXPECT_EQ(0, mag[2]);
//
//     field = gen_2d_heis_afm(0, 0, -1);
//     mag = hamil.calc_M(&field);
//     EXPECT_EQ(0, mag[0]);
//     EXPECT_EQ(0, mag[1]);
//     EXPECT_EQ(0, mag[2]);
// }
//
// TEST(Heis_model, 2d_submag)
// {
//     field_2d_h field = gen_2d_heis_fm(1, 0, 0);
//     ham_heis hamil(0, 1);
//     hamil.init_dim(&field);
//     std::vector<double> mag = hamil.calc_subM(&field, 1);
//     EXPECT_EQ(50, mag[0]);
//     EXPECT_EQ(0, mag[1]);
//     EXPECT_EQ(0, mag[2]);
//
//     field = gen_2d_heis_fm(-1, 0, 0);
//     mag = hamil.calc_subM(&field, 1);
//     EXPECT_EQ(-50, mag[0]);
//     EXPECT_EQ(0, mag[1]);
//     EXPECT_EQ(0, mag[2]);
//
//     field = gen_2d_heis_fm(0, 1, 0);
//     mag = hamil.calc_subM(&field, 1);
//     EXPECT_EQ(0, mag[0]);
//     EXPECT_EQ(50, mag[1]);
//     EXPECT_EQ(0, mag[2]);
//
//     field = gen_2d_heis_fm(0, -1, 0);
//     mag = hamil.calc_subM(&field, 1);
//     EXPECT_EQ(0, mag[0]);
//     EXPECT_EQ(-50, mag[1]);
//     EXPECT_EQ(0, mag[2]);
//
//     field = gen_2d_heis_fm(0, 0, 1);
//     mag = hamil.calc_subM(&field, 1);
//     EXPECT_EQ(0, mag[0]);
//     EXPECT_EQ(0, mag[1]);
//     EXPECT_EQ(50, mag[2]);
//
//     field = gen_2d_heis_fm(0, 0, -1);
//     mag = hamil.calc_subM(&field, 1);
//     EXPECT_EQ(0, mag[0]);
//     EXPECT_EQ(0, mag[1]);
//     EXPECT_EQ(-50, mag[2]);
//
//     field = gen_2d_heis_afm(1, 0, 0);
//     mag = hamil.calc_subM(&field, 1);
//     EXPECT_EQ(50, mag[0]);
//     EXPECT_EQ(0, mag[1]);
//     EXPECT_EQ(0, mag[2]);
//
//     field = gen_2d_heis_afm(-1, 0, 0);
//     mag = hamil.calc_subM(&field, 1);
//     EXPECT_EQ(-50, mag[0]);
//     EXPECT_EQ(0, mag[1]);
//     EXPECT_EQ(0, mag[2]);
//
//     field = gen_2d_heis_afm(0, 1, 0);
//     mag = hamil.calc_subM(&field, 1);
//     EXPECT_EQ(0, mag[0]);
//     EXPECT_EQ(50, mag[1]);
//     EXPECT_EQ(0, mag[2]);
//
//     field = gen_2d_heis_afm(0, -1, 0);
//     mag = hamil.calc_subM(&field, 1);
//     EXPECT_EQ(0, mag[0]);
//     EXPECT_EQ(-50, mag[1]);
//     EXPECT_EQ(0, mag[2]);
//
//     field = gen_2d_heis_afm(0, 0, 1);
//     mag = hamil.calc_subM(&field, 1);
//     EXPECT_EQ(0, mag[0]);
//     EXPECT_EQ(0, mag[1]);
//     EXPECT_EQ(50, mag[2]);
//
//     field = gen_2d_heis_afm(0, 0, -1);
//     mag = hamil.calc_subM(&field, 1);
//     EXPECT_EQ(0, mag[0]);
//     EXPECT_EQ(0, mag[1]);
//     EXPECT_EQ(-50, mag[2]);
// }
//
// TEST(Heis_model, 2d_dE_consist)
// {
//     field_2d_h field = gen_2d_heis_fm(1, 0, 0);
//     ham_heis hamil(0, 1);
//     hamil.init_dim(&field);
//     std::vector<int> pos(2);
//     double old_E = hamil.calc_E(&field);
//     for(int i = 0; i < 1000; i++)
//     {
//         pos[0] = int(st_rand_double.gen()*10 + 1);
//         pos[1] = int(st_rand_double.gen()*10 + 1);
//         double dE = hamil.dE(&field, pos);
//         field.change_to_test(pos, &hamil);
//         double new_E = hamil.calc_E(&field);
//         EXPECT_FLOAT_EQ(old_E + dE, new_E);
//         old_E = new_E;
//     }
// }
//
// TEST(Heis_model, 2d_top_charge)
// {
//     field_2d_h field = gen_2d_skyrm();
//     ham_heis hamil(0, 1);
//     hamil.init_dim(&field);
//
//     std::vector<double> tcs = hamil.calc_top_charge(&field);
//
//     EXPECT_FLOAT_EQ(tcs[0], -1);
// }

///////////////////////////////////////////////////////
// Heisenberg model tests - 3D
///////////////////////////////////////////////////////

TEST(Heis_model, 3d_energy_zero_field)
{
    EXPECT_DOUBLE_EQ(-2700, exchangeOnly.calc_E(heisFMField, zeroH));

    EXPECT_DOUBLE_EQ(2700, exchangeOnly.calc_E(heisAFMField, zeroH));
}

TEST(Heis_model, 3d_energy_ext_field)
{
    EXPECT_DOUBLE_EQ(-2800, exchangeOnly.calc_E(heisFMField, heisH));

    EXPECT_DOUBLE_EQ(2700, exchangeOnly.calc_E(heisAFMField, heisH));
}

TEST(Heis_model, 3d_mag)
{
    EXPECT_DOUBLE_EQ(0, exchangeOnly.calc_M(heisFMField)[0]);

    EXPECT_DOUBLE_EQ(0, exchangeOnly.calc_M(heisAFMField)[0]);

    EXPECT_DOUBLE_EQ(0, exchangeOnly.calc_M(heisFMField)[1]);

    EXPECT_DOUBLE_EQ(0, exchangeOnly.calc_M(heisAFMField)[1]);

    EXPECT_DOUBLE_EQ(1000, exchangeOnly.calc_M(heisFMField)[2]);

    EXPECT_DOUBLE_EQ(0, exchangeOnly.calc_M(heisAFMField)[2]);
}

TEST(Heis_model, 3d_submag)
{
    EXPECT_DOUBLE_EQ(0, exchangeOnly.calc_subM(heisFMField, 0)[0]);

    EXPECT_DOUBLE_EQ(0, exchangeOnly.calc_subM(heisAFMField, 0)[0]);

    EXPECT_DOUBLE_EQ(0, exchangeOnly.calc_subM(heisFMField, 0)[1]);

    EXPECT_DOUBLE_EQ(0, exchangeOnly.calc_subM(heisAFMField, 0)[1]);

    EXPECT_DOUBLE_EQ(500, exchangeOnly.calc_subM(heisFMField, 0)[2]);

    EXPECT_DOUBLE_EQ(500, exchangeOnly.calc_subM(heisAFMField, 0)[2]);
}

TEST(Heis_model, 3d_sub4mag)
{
    EXPECT_DOUBLE_EQ(0, exchangeOnly.calc_sub4M(heisFMField)[0]);

    EXPECT_DOUBLE_EQ(0, exchangeOnly.calc_sub4M(heisAFMField)[0]);

    EXPECT_DOUBLE_EQ(0, exchangeOnly.calc_sub4M(heisFMField)[1]);

    EXPECT_DOUBLE_EQ(0, exchangeOnly.calc_sub4M(heisAFMField)[1]);

    EXPECT_DOUBLE_EQ(63, exchangeOnly.calc_sub4M(heisFMField)[2]);

    EXPECT_DOUBLE_EQ(63, exchangeOnly.calc_sub4M(heisAFMField)[2]);
}

TEST(Heis_model, 3d_dE_consist)
{
    particle::field::field_type fieldCopy = heisFMField;
    int pos=0;
    double old_E = exchangeOnly.calc_E(fieldCopy, heisH);
    for(int i = 0; i < 1000; i++)
    {
        pos = int(st_rand_double.gen()*1000);
        fieldCopy.gen_rand();
        double dE = exchangeOnly.calc_dE(fieldCopy, pos, heisH);
        fieldCopy.set_rand(pos);
        double new_E = exchangeOnly.calc_E(fieldCopy, heisH);
        EXPECT_NEAR(old_E + dE, new_E, 1e-10);
        old_E = new_E;
    }
}

TEST(Heis_model, 3d_dE_consist_skyrm)
{
    particle::field::field_type fieldCopy = skyrmField;
    int pos=0;
    double old_E = exchangeOnly.calc_E(fieldCopy, heisH);
    for(int i = 0; i < 1000; i++)
    {
        pos = int(st_rand_double.gen()*1000);
        fieldCopy.gen_rand();
        double dE = exchangeOnly.calc_dE(fieldCopy, pos, heisH);
        fieldCopy.set_rand(pos);
        double new_E = exchangeOnly.calc_E(fieldCopy, heisH);
        EXPECT_NEAR(old_E + dE, new_E, 1e-10);
        old_E = new_E;
    }
}

TEST(Heis_model, 3d_top_charge)
{
    std::vector<double> tcs = exchangeOnly.calc_TC(heisFMField);

    for(int i = 0; i < tcs.size(); i++)
    {
        EXPECT_FLOAT_EQ(tcs[i], 0);
    }

    tcs = exchangeOnly.calc_TC(skyrmField);

    for(int i = 0; i < tcs.size(); i++)
    {
        EXPECT_NEAR(tcs[i], -1, 0.01);
    }
}

#endif
