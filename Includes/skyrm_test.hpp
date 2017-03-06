#ifndef _SKYRTEST
#define _SKYRTEST

#include "test_functions.hpp"

///////////////////////////////////////////////////////
// Skyrmion model tests - 3D
///////////////////////////////////////////////////////

TEST(Skyrmion, 3d_energy_zero_field)
{
    field_3d_h field = gen_3d_heis_fm(1, 0, 0);
    ham_skyrm hamil(0, 1, 0.75);
    hamil.init_dim(&field);
    // EXPECT_EQ(-2700, hamil.calc_E(&field));
    //
    // field = gen_3d_heis_fm(-1, 0, 0);
    // EXPECT_EQ(-2700, hamil.calc_E(&field));
    //
    // field = gen_3d_heis_fm(0, 1, 0);
    // EXPECT_EQ(-2700, hamil.calc_E(&field));
    //
    // field = gen_3d_heis_fm(0, -1, 0);
    // EXPECT_EQ(-2700, hamil.calc_E(&field));
    //
    // field = gen_3d_heis_fm(0, 0, 1);
    // EXPECT_EQ(-2700, hamil.calc_E(&field));
    //
    // field = gen_3d_heis_fm(0, 0, -1);
    // EXPECT_EQ(-2700, hamil.calc_E(&field));
    //
    // field = gen_3d_heis_afm(1, 0, 0);
    // EXPECT_EQ(2700, hamil.calc_E(&field));
    //
    // field = gen_3d_heis_afm(-1, 0, 0);
    // EXPECT_EQ(2700, hamil.calc_E(&field));
    //
    // field = gen_3d_heis_afm(0, 1, 0);
    // EXPECT_EQ(2700, hamil.calc_E(&field));
    //
    // field = gen_3d_heis_afm(0, -1, 0);
    // EXPECT_EQ(2700, hamil.calc_E(&field));
    //
    // field = gen_3d_heis_afm(0, 0, 1);
    // EXPECT_EQ(2700, hamil.calc_E(&field));
    //
    // field = gen_3d_heis_afm(0, 0, -1);
    // EXPECT_EQ(2700, hamil.calc_E(&field));

    vector<int> pos(3);
    for(int i = 1; i < 11; i++)
    {
        pos[0] = i;
        for(int j = 1; j < 11; j++)
        {
            pos[1] = j;
            for(int k = 1; k < 11; k++)
            {
                pos[2] = k;
                // cout << (i+j+k)%2 << " " << ((i+j+k)%2-1)*(-1) << endl;
                field.fill_val_h(pos, (i+j+k)%2, ((i+j+k)%2-1)*(-1), 0);
            }
        }
    }
    EXPECT_EQ(-900, hamil.calc_E(&field));
}

// TEST(Skyrmion, 3d_energy_ext_field)
// {
//     field_3d_h field = gen_3d_heis_fm(1, 0, 0);
//     ham_skyrm hamil(0.1, 1, 0.75);
//     hamil.init_dim(&field);
//     EXPECT_EQ(-2700, hamil.calc_E(&field));
//
//     field = gen_3d_heis_fm(-1, 0, 0);
//     EXPECT_EQ(-2700, hamil.calc_E(&field));
//
//     field = gen_3d_heis_fm(0, 1, 0);
//     EXPECT_EQ(-2700, hamil.calc_E(&field));
//
//     field = gen_3d_heis_fm(0, -1, 0);
//     EXPECT_EQ(-2700, hamil.calc_E(&field));
//
//     field = gen_3d_heis_fm(0, 0, 1);
//     EXPECT_EQ(-2800, hamil.calc_E(&field));
//
//     field = gen_3d_heis_fm(0, 0, -1);
//     EXPECT_EQ(-2600, hamil.calc_E(&field));
//
//     field = gen_3d_heis_afm(1, 0, 0);
//     EXPECT_EQ(2700, hamil.calc_E(&field));
//
//     field = gen_3d_heis_afm(-1, 0, 0);
//     EXPECT_EQ(2700, hamil.calc_E(&field));
//
//     field = gen_3d_heis_afm(0, 1, 0);
//     EXPECT_EQ(2700, hamil.calc_E(&field));
//
//     field = gen_3d_heis_afm(0, -1, 0);
//     EXPECT_EQ(2700, hamil.calc_E(&field));
//
//     field = gen_3d_heis_afm(0, 0, 1);
//     EXPECT_EQ(2700, hamil.calc_E(&field));
//
//     field = gen_3d_heis_afm(0, 0, -1);
//     EXPECT_EQ(2700, hamil.calc_E(&field));
// }

TEST(Skyrmion, 3d_mag)
{
    field_3d_h field = gen_3d_heis_fm(1, 0, 0);
    ham_skyrm hamil(0, 1, 0.75);
    hamil.init_dim(&field);
    vector<double> mag = hamil.calc_M(&field);
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

TEST(Skyrmion, 3d_submag)
{
    field_3d_h field = gen_3d_heis_fm(1, 0, 0);
    ham_skyrm hamil(0, 1, 0.75);
    hamil.init_dim(&field);
    vector<double> mag = hamil.calc_subM(&field, 1);
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

// TEST(Skyrmion, 3d_dE_consist)
// {
//     field_3d_h field = gen_3d_heis_fm(1, 0, 0);
//     ham_skyrm hamil(0, 1, 0.75);
//     hamil.init_dim(&field);
//     vector<int> pos(3);
//     double old_E = hamil.calc_E(&field);
//     for(int i = 0; i < 1000; i++)
//     {
//         pos[0] = int(st_rand_double.gen()*10 + 1);
//         pos[1] = int(st_rand_double.gen()*10 + 1);
//         pos[2] = int(st_rand_double.gen()*10 + 1);
//         double dE = hamil.dE(&field, pos);
//         field.change_to_test(pos, &hamil);
//         double new_E = hamil.calc_E(&field);
//         EXPECT_NEAR(old_E + dE, new_E, abs(new_E*1e-10));
//         old_E = new_E;
//     }
// }

#endif
