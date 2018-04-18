#ifndef _ISINGTEST
#define _ISINGTEST

#include "test_functions.hpp"
#include "../includes/td_funcs.hpp"

///////////////////////////////////////////////////////
// Ising model tests - 2D
///////////////////////////////////////////////////////

// TEST(Ising_model, 2d_ferromagnetic_energy_zero_field)
// {
//     field_2d_i field = gen_2d_ising_fm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(-180, hamil.calc_E(&field));
// }
//
// TEST(Ising_model, 2d_ferromagnetic_energy_ext_field)
// {
//     field_2d_i field = gen_2d_ising_fm();
//     ham_ising hamil(0.1, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(-190, hamil.calc_E(&field));
// }
//
// TEST(Ising_model, 2d_ferromagnetic_mag)
// {
//     field_2d_i field = gen_2d_ising_fm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(100, hamil.calc_M(&field)[0]);
// }
//
// TEST(Ising_model, 2d_ferromagnetic_submag)
// {
//     field_2d_i field = gen_2d_ising_fm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(50, hamil.calc_subM(&field, 1)[0]);
// }
//
// TEST(Ising_model, 2d_ferromagnetic_dE)
// {
//     field_2d_i field = gen_2d_ising_fm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     std::vector<int> pos(2,5);
//     EXPECT_EQ(8, hamil.dE(&field, pos));
// }
//
// TEST(Ising_model, 2d_antiferromagnetic_energy_zero_field)
// {
//     field_2d_i field = gen_2d_ising_afm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(180, hamil.calc_E(&field));
// }
//
// TEST(Ising_model, 2d_antiferromagnetic_energy_ext_field)
// {
//     field_2d_i field = gen_2d_ising_afm();
//     ham_ising hamil(0.3, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(180, hamil.calc_E(&field));
// }
//
// TEST(Ising_model, 2d_antiferromagnetic_mag)
// {
//     field_2d_i field = gen_2d_ising_afm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(0, hamil.calc_M(&field)[0]);
// }
//
// TEST(Ising_model, 2d_antiferromagnetic_submag)
// {
//     field_2d_i field = gen_2d_ising_afm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(50, hamil.calc_subM(&field, 1)[0]);
// }
//
// TEST(Ising_model, 2d_antiferromagnetic_dE)
// {
//     field_2d_i field = gen_2d_ising_afm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     std::vector<int> pos(2,5);
//     EXPECT_EQ(-8, hamil.dE(&field, pos));
// }
//
// TEST(Ising_model, 2d_dE_consist)
// {
//     field_2d_i field = gen_2d_ising_fm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     std::vector<int> pos(2);
//     int old_E = hamil.calc_E(&field);
//     for(int i = 0; i < 1000; i++)
//     {
//         pos[0] = int(st_rand_double.gen()*10 + 1);
//         pos[1] = int(st_rand_double.gen()*10 + 1);
//         int dE = hamil.dE(&field, pos);
//         field.change_to_test(pos, &hamil);
//         int new_E = hamil.calc_E(&field);
//         EXPECT_EQ(old_E + dE, new_E);
//         old_E = new_E;
//     }
// }
//
// ///////////////////////////////////////////////////////
// // Ising model tests - 3D
// ///////////////////////////////////////////////////////

TEST(Ising_model, 3d_ferromagnetic_energy_zero_field)
{
    std::valarray<double> H = {0, 0, 0, 0};
    particle::field::field_type field = gen_fm(3, true, 1, 0);
    EXPECT_DOUBLE_EQ(-2700, particle::funcs::calc_E(field, H));
}

TEST(Ising_model, 3d_ferromagnetic_energy_ext_field)
{
    std::valarray<double> H = {0.1, 0, 0, 0};
    particle::field::field_type field = gen_fm(3, true, 1, 0);
    EXPECT_DOUBLE_EQ(-2800, particle::funcs::calc_E(field, H));
}

TEST(Ising_model, 3d_ferromagnetic_mag)
{
    particle::field::field_type field = gen_fm(3, true, 1, 0);
    EXPECT_DOUBLE_EQ(1000, particle::funcs::calc_M(field)[0]);
}

TEST(Ising_model, 3d_ferromagnetic_submag)
{
    particle::field::field_type field = gen_fm(3, true, 1, 0);
    EXPECT_DOUBLE_EQ(500, particle::funcs::calc_subM(field, 0)[0]);
}

TEST(Ising_model, 3d_ferromagnetic_dE)
{
    std::valarray<double> H = {0, 0, 0, 0};
    particle::field::field_type field = gen_fm(3, true, 1, 0);
    EXPECT_DOUBLE_EQ(12, particle::funcs::calc_dE(field, 555, H));
}

// TEST(Ising_model, 3d_antiferromagnetic_energy_zero_field)
// {
//     field_3d_i field = gen_3d_ising_afm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(2700, hamil.calc_E(&field));
// }
//
// TEST(Ising_model, 3d_antiferromagnetic_energy_ext_field)
// {
//     field_3d_i field = gen_3d_ising_afm();
//     ham_ising hamil(0.1, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(2700, hamil.calc_E(&field));
// }
//
// TEST(Ising_model, 3d_antiferromagnetic_mag)
// {
//     field_3d_i field = gen_3d_ising_afm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(0, hamil.calc_M(&field)[0]);
// }
//
// TEST(Ising_model, 3d_antiferromagnetic_submag)
// {
//     field_3d_i field = gen_3d_ising_fm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     EXPECT_EQ(500, hamil.calc_subM(&field, 1)[0]);
// }
//
// TEST(Ising_model, 3d_antiferromagnetic_dE)
// {
//     field_3d_i field = gen_3d_ising_fm();
//     ham_ising hamil(0, 1);
//     hamil.init_dim(&field);
//     std::vector<int> pos(3,5);
//     EXPECT_EQ(12, hamil.dE(&field, pos));
// }

TEST(Ising_model, 3d_dE_consist)
{
    std::valarray<double> H = {0, 0, 0, 0};
    particle::field::field_type field = gen_fm(3, true, 1, 0);
    int pos=0;
    double old_E = particle::funcs::calc_E(field, H);
    for(int i = 0; i < 1000; i++)
    {
        pos = int(st_rand_double.gen()*1000);
        double dE = particle::funcs::calc_dE(field, pos, H);
        field.set_rand(pos);
        // double new_E = particle::funcs::calc_E(field, H);
        // EXPECT_DOUBLE_EQ(old_E + dE, new_E);
        // old_E = new_E;
    }
}

#endif
