#ifndef _HEISTEST
#define _HEISTEST

#include "test_functions.hpp"

///////////////////////////////////////////////////////
// Heisenberg model tests - 2D
///////////////////////////////////////////////////////

TEST(Heis_model, 2d_energy_zero_field)
{
    field_2d_h field = gen_2d_heis_fm(1, 0, 0);
    ham_heis hamil(0, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(-180, hamil.calc_E(&field));

    field = gen_2d_heis_fm(-1, 0, 0);
    EXPECT_EQ(-180, hamil.calc_E(&field));

    field = gen_2d_heis_fm(0, 1, 0);
    EXPECT_EQ(-180, hamil.calc_E(&field));

    field = gen_2d_heis_fm(0, -1, 0);
    EXPECT_EQ(-180, hamil.calc_E(&field));

    field = gen_2d_heis_fm(0, 0, 1);
    EXPECT_EQ(-180, hamil.calc_E(&field));

    field = gen_2d_heis_fm(0, 0, -1);
    EXPECT_EQ(-180, hamil.calc_E(&field));

    field = gen_2d_heis_afm(1, 0, 0);
    EXPECT_EQ(180, hamil.calc_E(&field));

    field = gen_2d_heis_afm(-1, 0, 0);
    EXPECT_EQ(180, hamil.calc_E(&field));

    field = gen_2d_heis_afm(0, 1, 0);
    EXPECT_EQ(180, hamil.calc_E(&field));

    field = gen_2d_heis_afm(0, -1, 0);
    EXPECT_EQ(180, hamil.calc_E(&field));

    field = gen_2d_heis_afm(0, 0, 1);
    EXPECT_EQ(180, hamil.calc_E(&field));

    field = gen_2d_heis_afm(0, 0, -1);
    EXPECT_EQ(180, hamil.calc_E(&field));
}

TEST(Heis_model, 2d_energy_ext_field)
{
    field_2d_h field = gen_2d_heis_fm(1, 0, 0);
    ham_heis hamil(0.1, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(-180, hamil.calc_E(&field));

    field = gen_2d_heis_fm(-1, 0, 0);
    EXPECT_EQ(-180, hamil.calc_E(&field));

    field = gen_2d_heis_fm(0, 1, 0);
    EXPECT_EQ(-180, hamil.calc_E(&field));

    field = gen_2d_heis_fm(0, -1, 0);
    EXPECT_EQ(-180, hamil.calc_E(&field));

    field = gen_2d_heis_fm(0, 0, 1);
    EXPECT_EQ(-190, hamil.calc_E(&field));

    field = gen_2d_heis_fm(0, 0, -1);
    EXPECT_EQ(-170, hamil.calc_E(&field));

    field = gen_2d_heis_afm(1, 0, 0);
    EXPECT_EQ(180, hamil.calc_E(&field));

    field = gen_2d_heis_afm(-1, 0, 0);
    EXPECT_EQ(180, hamil.calc_E(&field));

    field = gen_2d_heis_afm(0, 1, 0);
    EXPECT_EQ(180, hamil.calc_E(&field));

    field = gen_2d_heis_afm(0, -1, 0);
    EXPECT_EQ(180, hamil.calc_E(&field));

    field = gen_2d_heis_afm(0, 0, 1);
    EXPECT_EQ(180, hamil.calc_E(&field));

    field = gen_2d_heis_afm(0, 0, -1);
    EXPECT_EQ(180, hamil.calc_E(&field));
}

TEST(Heis_model, 2d_mag)
{
    field_2d_h field = gen_2d_heis_fm(1, 0, 0);
    ham_heis hamil(0, 1);
    hamil.init_dim(&field);
    std::vector<double> mag = hamil.calc_M(&field);
    EXPECT_EQ(100, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_fm(-1, 0, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(-100, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_fm(0, 1, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(100, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_fm(0, -1, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(-100, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_fm(0, 0, 1);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(100, mag[2]);

    field = gen_2d_heis_fm(0, 0, -1);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(-100, mag[2]);

    field = gen_2d_heis_afm(1, 0, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_afm(-1, 0, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_afm(0, 1, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_afm(0, -1, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_afm(0, 0, 1);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_afm(0, 0, -1);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);
}

TEST(Heis_model, 2d_submag)
{
    field_2d_h field = gen_2d_heis_fm(1, 0, 0);
    ham_heis hamil(0, 1);
    hamil.init_dim(&field);
    std::vector<double> mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(50, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_fm(-1, 0, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(-50, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_fm(0, 1, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(50, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_fm(0, -1, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(-50, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_fm(0, 0, 1);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(50, mag[2]);

    field = gen_2d_heis_fm(0, 0, -1);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(-50, mag[2]);

    field = gen_2d_heis_afm(1, 0, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(50, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_afm(-1, 0, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(-50, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_afm(0, 1, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(50, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_afm(0, -1, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(-50, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_2d_heis_afm(0, 0, 1);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(50, mag[2]);

    field = gen_2d_heis_afm(0, 0, -1);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(-50, mag[2]);
}

TEST(Heis_model, 2d_dE_consist)
{
    field_2d_h field = gen_2d_heis_fm(1, 0, 0);
    ham_heis hamil(0, 1);
    hamil.init_dim(&field);
    std::vector<int> pos(2);
    double old_E = hamil.calc_E(&field);
    for(int i = 0; i < 1000; i++)
    {
        pos[0] = int(st_rand_double.gen()*10 + 1);
        pos[1] = int(st_rand_double.gen()*10 + 1);
        double dE = hamil.dE(&field, pos);
        field.change_to_test(pos, &hamil);
        double new_E = hamil.calc_E(&field);
        EXPECT_FLOAT_EQ(old_E + dE, new_E);
        old_E = new_E;
    }
}

TEST(Heis_model, 2d_top_charge)
{
    field_2d_h field = gen_2d_skyrm();
    ham_heis hamil(0, 1);
    hamil.init_dim(&field);

    std::vector<double> tcs = hamil.calc_top_charge(&field);

    EXPECT_FLOAT_EQ(tcs[0], -1);
}

///////////////////////////////////////////////////////
// Heisenberg model tests - 3D
///////////////////////////////////////////////////////

TEST(Heis_model, 3d_energy_zero_field)
{
    field_3d_h field = gen_3d_heis_fm(1, 0, 0);
    ham_heis hamil(0, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(-2700, hamil.calc_E(&field));

    field = gen_3d_heis_fm(-1, 0, 0);
    EXPECT_EQ(-2700, hamil.calc_E(&field));

    field = gen_3d_heis_fm(0, 1, 0);
    EXPECT_EQ(-2700, hamil.calc_E(&field));

    field = gen_3d_heis_fm(0, -1, 0);
    EXPECT_EQ(-2700, hamil.calc_E(&field));

    field = gen_3d_heis_fm(0, 0, 1);
    EXPECT_EQ(-2700, hamil.calc_E(&field));

    field = gen_3d_heis_fm(0, 0, -1);
    EXPECT_EQ(-2700, hamil.calc_E(&field));

    field = gen_3d_heis_afm(1, 0, 0);
    EXPECT_EQ(2700, hamil.calc_E(&field));

    field = gen_3d_heis_afm(-1, 0, 0);
    EXPECT_EQ(2700, hamil.calc_E(&field));

    field = gen_3d_heis_afm(0, 1, 0);
    EXPECT_EQ(2700, hamil.calc_E(&field));

    field = gen_3d_heis_afm(0, -1, 0);
    EXPECT_EQ(2700, hamil.calc_E(&field));

    field = gen_3d_heis_afm(0, 0, 1);
    EXPECT_EQ(2700, hamil.calc_E(&field));

    field = gen_3d_heis_afm(0, 0, -1);
    EXPECT_EQ(2700, hamil.calc_E(&field));
}

TEST(Heis_model, 3d_energy_ext_field)
{
    field_3d_h field = gen_3d_heis_fm(1, 0, 0);
    ham_heis hamil(0.1, 1);
    hamil.init_dim(&field);
    EXPECT_EQ(-2700, hamil.calc_E(&field));

    field = gen_3d_heis_fm(-1, 0, 0);
    EXPECT_EQ(-2700, hamil.calc_E(&field));

    field = gen_3d_heis_fm(0, 1, 0);
    EXPECT_EQ(-2700, hamil.calc_E(&field));

    field = gen_3d_heis_fm(0, -1, 0);
    EXPECT_EQ(-2700, hamil.calc_E(&field));

    field = gen_3d_heis_fm(0, 0, 1);
    EXPECT_EQ(-2800, hamil.calc_E(&field));

    field = gen_3d_heis_fm(0, 0, -1);
    EXPECT_EQ(-2600, hamil.calc_E(&field));

    field = gen_3d_heis_afm(1, 0, 0);
    EXPECT_EQ(2700, hamil.calc_E(&field));

    field = gen_3d_heis_afm(-1, 0, 0);
    EXPECT_EQ(2700, hamil.calc_E(&field));

    field = gen_3d_heis_afm(0, 1, 0);
    EXPECT_EQ(2700, hamil.calc_E(&field));

    field = gen_3d_heis_afm(0, -1, 0);
    EXPECT_EQ(2700, hamil.calc_E(&field));

    field = gen_3d_heis_afm(0, 0, 1);
    EXPECT_EQ(2700, hamil.calc_E(&field));

    field = gen_3d_heis_afm(0, 0, -1);
    EXPECT_EQ(2700, hamil.calc_E(&field));
}

TEST(Heis_model, 3d_mag)
{
    field_3d_h field = gen_3d_heis_fm(1, 0, 0);
    ham_heis hamil(0, 1);
    hamil.init_dim(&field);
    std::vector<double> mag = hamil.calc_M(&field);
    EXPECT_EQ(1000, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_fm(-1, 0, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(-1000, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_fm(0, 1, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(1000, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_fm(0, -1, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(-1000, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_fm(0, 0, 1);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(1000, mag[2]);

    field = gen_3d_heis_fm(0, 0, -1);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(-1000, mag[2]);

    field = gen_3d_heis_afm(1, 0, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_afm(-1, 0, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_afm(0, 1, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_afm(0, -1, 0);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_afm(0, 0, 1);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_afm(0, 0, -1);
    mag = hamil.calc_M(&field);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);
}

TEST(Heis_model, 3d_submag)
{
    field_3d_h field = gen_3d_heis_fm(1, 0, 0);
    ham_heis hamil(0, 1);
    hamil.init_dim(&field);
    std::vector<double> mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(500, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_fm(-1, 0, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(-500, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_fm(0, 1, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(500, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_fm(0, -1, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(-500, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_fm(0, 0, 1);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(500, mag[2]);

    field = gen_3d_heis_fm(0, 0, -1);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(-500, mag[2]);

    field = gen_3d_heis_afm(1, 0, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(500, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_afm(-1, 0, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(-500, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_afm(0, 1, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(500, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_afm(0, -1, 0);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(-500, mag[1]);
    EXPECT_EQ(0, mag[2]);

    field = gen_3d_heis_afm(0, 0, 1);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(500, mag[2]);

    field = gen_3d_heis_afm(0, 0, -1);
    mag = hamil.calc_subM(&field, 1);
    EXPECT_EQ(0, mag[0]);
    EXPECT_EQ(0, mag[1]);
    EXPECT_EQ(-500, mag[2]);
}

TEST(Heis_model, 3d_dE_consist)
{
    field_3d_h field = gen_3d_heis_fm(1, 0, 0);
    ham_heis hamil(0, 1);
    hamil.init_dim(&field);
    std::vector<int> pos(3);
    double old_E = hamil.calc_E(&field);
    for(int i = 0; i < 1000; i++)
    {
        pos[0] = int(st_rand_double.gen()*10 + 1);
        pos[1] = int(st_rand_double.gen()*10 + 1);
        pos[2] = int(st_rand_double.gen()*10 + 1);
        double dE = hamil.dE(&field, pos);
        field.change_to_test(pos, &hamil);
        double new_E = hamil.calc_E(&field);
        EXPECT_FLOAT_EQ(old_E + dE, new_E);
        old_E = new_E;
    }
}

TEST(Heis_model, 3d_top_charge)
{
    field_3d_h field = gen_3d_skyrm();
    ham_heis hamil(0, 1);
    hamil.init_dim(&field);

    std::vector<double> tcs = hamil.calc_top_charge(&field);

    for(int i = 0; i < 50; i++)
    {
        EXPECT_FLOAT_EQ(tcs[i], -1);
    }
}

#endif
