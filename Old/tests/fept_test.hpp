#ifndef _FEPTTEST
#define _FEPTTEST

#include "test_functions.hpp"

///////////////////////////////////////////////////////
// FePt model tests - 3D
///////////////////////////////////////////////////////

TEST(FePt, 3d_energy_zero_field)
{
    field_3d_h field = gen_3d_heis_fm(1, 0, 0);
    ham_FePt hamil;
    hamil.init_dim(&field);
    EXPECT_FLOAT_EQ(-94845.01392843794, hamil.calc_E(&field));

    field = gen_3d_heis_fm(0, 1, 0);
    EXPECT_FLOAT_EQ(-94845.01392843794, hamil.calc_E(&field));

    field = gen_3d_heis_fm(0, 0, 1);
    EXPECT_FLOAT_EQ(-96045.87686933874, hamil.calc_E(&field));

    field = gen_3d_heis_fm(-1, 0, 0);
    EXPECT_FLOAT_EQ(-94845.01392843794, hamil.calc_E(&field));

    field = gen_3d_heis_fm(0, -1, 0);
    EXPECT_FLOAT_EQ(-94845.01392843794, hamil.calc_E(&field));

    field = gen_3d_heis_fm(0, 0, -1);
    EXPECT_FLOAT_EQ(-96045.87686933874, hamil.calc_E(&field));

    // Anti
    field = gen_3d_heis_afm(1, 0, 0);
    EXPECT_FLOAT_EQ(14349.191660944736, hamil.calc_E(&field));

    field = gen_3d_heis_afm(0, 1, 0);
    EXPECT_FLOAT_EQ(14349.191660944736, hamil.calc_E(&field));

    field = gen_3d_heis_afm(0, 0, 1);
    EXPECT_FLOAT_EQ(14481.605270298143, hamil.calc_E(&field));

    field = gen_3d_heis_afm(-1, 0, 0);
    EXPECT_FLOAT_EQ(14349.191660944736, hamil.calc_E(&field));

    field = gen_3d_heis_afm(0, -1, 0);
    EXPECT_FLOAT_EQ(14349.191660944736, hamil.calc_E(&field));

    field = gen_3d_heis_afm(0, 0, -1);
    EXPECT_FLOAT_EQ(14481.605270298143, hamil.calc_E(&field));
}

TEST(FePt, 3d_energy_ext_field)
{
    field_3d_h field = gen_3d_heis_fm(1, 0, 0);
    ham_FePt hamil(1);
    hamil.init_dim(&field);
    EXPECT_FLOAT_EQ(-94845.01392843794, hamil.calc_E(&field));

    field = gen_3d_heis_fm(0, 1, 0);
    EXPECT_FLOAT_EQ(-94845.01392843794, hamil.calc_E(&field));

    field = gen_3d_heis_fm(0, 0, 1);
    EXPECT_FLOAT_EQ(-97045.87686933874, hamil.calc_E(&field));

    field = gen_3d_heis_fm(-1, 0, 0);
    EXPECT_FLOAT_EQ(-94845.01392843794, hamil.calc_E(&field));

    field = gen_3d_heis_fm(0, -1, 0);
    EXPECT_FLOAT_EQ(-94845.01392843794, hamil.calc_E(&field));

    field = gen_3d_heis_fm(0, 0, -1);
    EXPECT_FLOAT_EQ(-95045.87686933874, hamil.calc_E(&field));

    // Anti
    field = gen_3d_heis_afm(1, 0, 0);
    EXPECT_FLOAT_EQ(14349.191660944736, hamil.calc_E(&field));

    field = gen_3d_heis_afm(0, 1, 0);
    EXPECT_FLOAT_EQ(14349.191660944736, hamil.calc_E(&field));

    field = gen_3d_heis_afm(0, 0, 1);
    EXPECT_FLOAT_EQ(14481.605270298143, hamil.calc_E(&field));

    field = gen_3d_heis_afm(-1, 0, 0);
    EXPECT_FLOAT_EQ(14349.191660944736, hamil.calc_E(&field));

    field = gen_3d_heis_afm(0, -1, 0);
    EXPECT_FLOAT_EQ(14349.191660944736, hamil.calc_E(&field));

    field = gen_3d_heis_afm(0, 0, -1);
    EXPECT_FLOAT_EQ(14481.605270298143, hamil.calc_E(&field));
}

TEST(FePt, 3d_mag)
{
    field_3d_h field = gen_3d_heis_fm(1, 0, 0);
    ham_FePt hamil;
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

TEST(FePt, 3d_submag)
{
    field_3d_h field = gen_3d_heis_fm(1, 0, 0);
    ham_FePt hamil;
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

TEST(FePt, 3d_dE_consist)
{
    field_3d_h field = gen_3d_heis_fm(1, 0, 0);
    ham_FePt hamil(4);
    hamil.init_dim(&field);
    std::vector<int> pos(3);
    double old_E = hamil.calc_E(&field);
    for(int i = 0; i < 100; i++)
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

#endif
