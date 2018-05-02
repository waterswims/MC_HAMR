#include "../includes/output.hpp"
#include "../includes/array_alloc.hpp"
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <hdf5.h>

std::string check_h5_file(const double J,
                           const double size,
                           const double am,
                           const double asd,
                           const double k,
                           const double K,
                           const char shape,
                           const char hamil,
                           const bool perio,
                           const int prot,
                           const bool distrib,
                           const int N_temps,
                           const int N_fields,
                           const int N_samps,
                           const int N_latts,
                           const float* Ts,
                           const float* Hs,
                           const int v1_size,
                           const int tc_size,
                           bool** checkp,
                           bool &file_exists)
{
   std::stringstream prestream, datstream;
   std::string datname;
   std::string prefix;
   prestream << "Output/" << "J_" << J;
   if (distrib)
   {
       prestream << "-am_" << am << "-asd_" << asd;
   }
   else
   {
       prestream << "-s_" << size;
   }
   prestream << "-k_" << k;
   if (hamil == 's' || hamil == 'S')
   {
       prestream << "-K_" << K;
   }
   prestream << "-sh_" << shape << "-ha_" << hamil << "-per_" << int(perio);
   prestream << "-pr_" << prot << ".h5";
   prestream >> prefix;

   // Check the file exists
   std::ifstream fin(prefix);
   if (!fin)
   {
    //    std::cout << "Here" << std::endl;
    //    exit(101);
       create_h5_file(prefix, hamil, distrib, N_temps, N_fields, N_samps,
                      N_latts, Ts, Hs, v1_size, tc_size);
       for(int i=0; i<v1_size; i++)
       {
           for(int j=0; j<N_latts; j++)
           {
               checkp[i][j] = false;
           }
       }
   }
   else
   {
       fin.close();
       // Open file
       hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
       H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
       hid_t f_id = H5Fopen(prefix.c_str(), H5F_ACC_RDONLY, plist_id);
       H5Pclose(plist_id);
       // Read file
       hid_t dset_id = H5Dopen1(f_id, "Complete");
       int* input = alloc_1darr<int>(v1_size*N_latts);
       H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, input);
       int k=0;
       for(int i=0; i<v1_size; i++)
       {
           for(int j=0; j<N_latts; j++)
           {
               checkp[i][j] = input[k];
               k++;
           }
       }
       // Close file
       dealloc_1darr<int>(input);
       H5Dclose(dset_id);
       H5Fclose(f_id);
       file_exists = true;
   }
   return prefix;
}

void create_h5_file(std::string prefix,
                   const char hamil,
                   const bool distrib,
                   const int N_temps,
                   const int N_fields,
                   const int N_samps,
                   const int N_latts,
                   const float* Ts,
                   const float* Hs,
                   const int v1_size,
                   const int tc_size)
{
    // create file
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t file_id;
    file_id = H5Fcreate(prefix.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    // create fulldat group
    hid_t g_id;
    g_id = H5Gcreate2(file_id, "/Fulldat", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(g_id);

    if (!distrib)
    {
        // Put in latt print group
        g_id = H5Gcreate2(file_id, "/Av_Latt", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(g_id);

        g_id = H5Gcreate2(file_id, "/Sing_Latt", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(g_id);
    }

    // create datasets
    hsize_t pd_dims[3] = {N_fields, N_temps, N_samps*N_latts};
    hid_t dspace_id = H5Screate_simple(3, pd_dims, NULL);
    hid_t dset_id = H5Dcreate(file_id, "/Fulldat/mags", H5T_NATIVE_FLOAT,
        dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dclose(dset_id);
    dset_id = H5Dcreate(file_id, "/Fulldat/energies", H5T_NATIVE_FLOAT,
        dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dclose(dset_id);
    dset_id = H5Dcreate(file_id, "/Fulldat/sub_mags", H5T_NATIVE_FLOAT,
        dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dclose(dset_id);
    if (hamil != 'i' && hamil != 'I')
    {
        dset_id = H5Dcreate(file_id, "/Fulldat/mag_xs", H5T_NATIVE_FLOAT,
            dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dclose(dset_id);
        dset_id = H5Dcreate(file_id, "/Fulldat/mag_ys", H5T_NATIVE_FLOAT,
            dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dclose(dset_id);
        dset_id = H5Dcreate(file_id, "/Fulldat/mag_zs", H5T_NATIVE_FLOAT,
            dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dclose(dset_id);
        dset_id = H5Dcreate(file_id, "/Fulldat/sub_mag_xs", H5T_NATIVE_FLOAT,
            dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dclose(dset_id);
        dset_id = H5Dcreate(file_id, "/Fulldat/sub_mag_ys", H5T_NATIVE_FLOAT,
            dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dclose(dset_id);
        dset_id = H5Dcreate(file_id, "/Fulldat/sub_mag_zs", H5T_NATIVE_FLOAT,
            dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dclose(dset_id);
    }
    H5Sclose(dspace_id);
    hsize_t num_dims[1] = {N_samps};
    dspace_id = H5Screate_simple(1, num_dims, NULL);
    dset_id = H5Dcreate(file_id, "/Fulldat/nums", H5T_NATIVE_FLOAT,
        dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dclose(dset_id);
    H5Sclose(dspace_id);

    num_dims[0] = N_temps;
    dspace_id = H5Screate_simple(1, num_dims, NULL);
    dset_id = H5Dcreate(file_id, "/Ts", H5T_NATIVE_FLOAT,
        dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dwrite (dset_id, H5T_NATIVE_FLOAT, dspace_id, dspace_id, plist_id, Ts);
    H5Dclose(dset_id);
    H5Sclose(dspace_id);
    H5Pclose(plist_id);

    num_dims[0] = N_fields;
    dspace_id = H5Screate_simple(1, num_dims, NULL);
    dset_id = H5Dcreate(file_id, "/Hs", H5T_NATIVE_FLOAT,
        dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dwrite (dset_id, H5T_NATIVE_FLOAT, dspace_id, dspace_id, plist_id, Hs);
    H5Dclose(dset_id);
    H5Sclose(dspace_id);
    H5Pclose(plist_id);

    hsize_t check_dims[2] = {v1_size, N_latts};
    dspace_id = H5Screate_simple(2, check_dims, NULL);
    dset_id = H5Dcreate(file_id, "/Complete", H5T_NATIVE_INT,
        dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dclose(dset_id);
    H5Sclose(dspace_id);

    if (hamil != 'i' && hamil != 'I' && !distrib)
    {
        hsize_t tc_dims[4] = {N_fields, N_temps, tc_size, N_samps*N_latts};
        dspace_id = H5Screate_simple(4, tc_dims, NULL);
        dset_id = H5Dcreate(file_id, "/Fulldat/top_chars", H5T_NATIVE_FLOAT,
            dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dclose(dset_id);
        H5Sclose(dspace_id);
    }

    H5Fclose(file_id);
}

void print_TD_h5(const float* magx,
             const float* magy,
             const float* magz,
             const float* mag,
             const float* ener,
             const float* smagx,
             const float* smagy,
             const float* smagz,
             const float* smag,
             float** tcs,
             const int N_samp,
             const int protocol,
             const int var1,
             const int var2,
             const int v2max,
             const int tc_size,
             const std::string prefix,
             const char hamil,
             const bool distrib,
             const int latt_num)
{
    // Get protocol position
    int H, T;
    switch(protocol)
    {
        case 1:
        H = var1;
        T = var2;
        break;
        case 2:
        case 4:
        H = var2;
        T = var1;
        break;
    }

    // Open the file
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t f_id = H5Fopen(prefix.c_str(), H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);

    // Open and write mag
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    hid_t dset_id = H5Dopen1(f_id, "/Fulldat/mags");
    hid_t slab_id = H5Dget_space(dset_id);
    hsize_t count[3] = {1, 1, N_samp};
    hsize_t offset[3] = {H, T, latt_num*N_samp};
    H5Sselect_hyperslab(slab_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    hid_t dspace_id = H5Screate_simple(3, count, NULL);
    H5Dwrite(dset_id, H5T_NATIVE_FLOAT, dspace_id, slab_id, plist_id, mag);
    H5Sclose(slab_id);
    H5Dclose(dset_id);

    // Other variables
    dset_id = H5Dopen1(f_id, "/Fulldat/energies");
    slab_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(slab_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    H5Dwrite(dset_id, H5T_NATIVE_FLOAT, dspace_id, slab_id, plist_id, ener);
    H5Sclose(slab_id);
    H5Dclose(dset_id);

    dset_id = H5Dopen1(f_id, "/Fulldat/sub_mags");
    slab_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(slab_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    H5Dwrite(dset_id, H5T_NATIVE_FLOAT, dspace_id, slab_id, plist_id, smag);
    H5Sclose(slab_id);
    H5Dclose(dset_id);

    if(hamil != 'i' || hamil != 'I')
    {
        dset_id = H5Dopen1(f_id, "/Fulldat/mag_xs");
        slab_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(slab_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Dwrite(dset_id, H5T_NATIVE_FLOAT, dspace_id, slab_id, plist_id,
                 magx);
        H5Sclose(slab_id);
        H5Dclose(dset_id);

        dset_id = H5Dopen1(f_id, "/Fulldat/mag_ys");
        slab_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(slab_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Dwrite(dset_id, H5T_NATIVE_FLOAT, dspace_id, slab_id, plist_id,
                 magy);
        H5Sclose(slab_id);
        H5Dclose(dset_id);

        dset_id = H5Dopen1(f_id, "/Fulldat/mag_zs");
        slab_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(slab_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Dwrite(dset_id, H5T_NATIVE_FLOAT, dspace_id, slab_id, plist_id,
                 magz);
        H5Sclose(slab_id);
        H5Dclose(dset_id);

        dset_id = H5Dopen1(f_id, "/Fulldat/sub_mag_xs");
        slab_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(slab_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Dwrite(dset_id, H5T_NATIVE_FLOAT, dspace_id, slab_id, plist_id,
                 smagx);
        H5Sclose(slab_id);
        H5Dclose(dset_id);

        dset_id = H5Dopen1(f_id, "/Fulldat/sub_mag_ys");
        slab_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(slab_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Dwrite(dset_id, H5T_NATIVE_FLOAT, dspace_id, slab_id, plist_id,
                 smagy);
        H5Sclose(slab_id);
        H5Dclose(dset_id);

        dset_id = H5Dopen1(f_id, "/Fulldat/sub_mag_zs");
        slab_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(slab_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Dwrite(dset_id, H5T_NATIVE_FLOAT, dspace_id, slab_id, plist_id,
                 smagz);
        H5Sclose(slab_id);
        H5Dclose(dset_id);
    }
    H5Sclose(dspace_id);

    if (hamil != 'i' && hamil != 'I' && !distrib)
    {
        hsize_t count3[4] = {1, 1, 1, N_samp};
        hsize_t offset3[4] = {H, T, 0, latt_num*N_samp};

        // Top charges
        dset_id = H5Dopen1(f_id, "/Fulldat/top_chars");
        for(int tc_idx = 0; tc_idx < tc_size; tc_idx++)
        {
            offset3[2] = tc_idx;
            slab_id = H5Dget_space(dset_id);
            dspace_id = H5Screate_simple(4, count3, NULL);
            H5Sselect_hyperslab(slab_id, H5S_SELECT_SET, offset3, NULL, count3,
                NULL);
            H5Dwrite(dset_id, H5T_NATIVE_FLOAT, dspace_id, slab_id, plist_id,
                tcs[tc_idx]);
            H5Sclose(slab_id);
            H5Sclose(dspace_id);
        }
        H5Dclose(dset_id);
    }

    // Print the checkpoint
    if(var2 == v2max - 1)
    {
        dset_id = H5Dopen1(f_id, "/Complete");
        slab_id = H5Dget_space(dset_id);
        hsize_t count2[2] = {1, 1};
        hsize_t offset2[2] = {var1, latt_num};
        H5Sselect_hyperslab(slab_id, H5S_SELECT_SET, offset2, NULL, count2, NULL);
        dspace_id = H5Screate_simple(2, count2, NULL);
        int cpoint_val = 1;
        H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_id, slab_id, plist_id,
                 &cpoint_val);
        H5Sclose(slab_id);
        H5Dclose(dset_id);
        H5Sclose(dspace_id);
    }

    // Close File
    H5Pclose(plist_id);
    H5Fclose(f_id);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        std::cout << "\rPrinting H = " << H << ", T = " << T;
        std::cout << ", lattice number " << latt_num << "              ";
    }
}
