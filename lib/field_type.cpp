#include "../includes/field_type.hpp"
#include "../includes/array_alloc.hpp"

#ifdef __INTEL_COMPILER
#include "../includes/mklrand.hpp"
#define IRANDTYPE mklrand::mkl_irand
#define DRANDTYPE mklrand::mkl_drand
#else
#include "../includes/stdrand.hpp"
#define IRANDTYPE stdrand::std_i_unirand
#define DRANDTYPE stdrand::std_d_unirand
#endif

#include <cmath>
#include <cstdlib>
#include <hdf5.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <mpi.h>

extern IRANDTYPE st_rand_int;
extern DRANDTYPE st_rand_double;

particle::field::field_type::field_type(bool ising_in,
    bool periodic_in,
    int d_in,
    int edgesize_in,
    double J_mod,
    double D_mod,
    std::string J_filename)
{
    ising = ising_in;
    periodic = periodic_in;
    d = d_in;
    edgesize = edgesize_in;

    this->set_default_spins();

    std::ifstream Jstream;
    Jstream.open(J_filename.c_str());
    int icurr;
    double dcurr;
    while(Jstream >> icurr)
    {
        std::valarray<int> new_loc = blankloc;
        std::valarray<double> new_dmi = upspin;

        new_loc[0] = icurr;
        for (int i = 1; i < d; i++)
        {
            Jstream >> icurr;
            new_loc[i] = icurr;
        }
        loc_diffs.push_back(new_loc);

        Jstream >> dcurr;
        J_diffs.push_back(dcurr);

        for (int i = 0; i < 3; i++)
        {
            Jstream >> dcurr;
            new_dmi[i] = dcurr;
        }
        D_vecs.push_back(new_dmi);
    }

    Jstream.close();
}

void particle::field::field_type::set_default_spins()
{
    blankloc.resize(d);
    int ssize = 4;
    upspin.resize(ssize);
    upspin = 0;
    downspin = upspin;
    testspin = upspin;

    if(ising)
    {
        upspin[0] = 1;
        downspin[0] = -1;
    }
    else
    {
        upspin[2] = 1;
        downspin[2] = -1;
    }

}

void particle::field::field_type::add_spin(std::valarray<int>& loc)
{
    spins.push_back(upspin);
    locs.push_back(loc);
}

void particle::field::field_type::set_neigh()
{
    std::vector<int> n_base;
    neighbours.resize(spins.size());
    neigh_choice.resize(spins.size());
    std::valarray<bool> check_arr;
    int nn = 0;
    for(unsigned int i = 0; i < spins.size(); i++)
    {
        neighbours[i] = n_base;
        neigh_choice[i] = n_base;
        for(unsigned int j = 0; j < spins.size(); j++)
        {
            if(i == j) {continue;}
            for(unsigned int k = 0; k < loc_diffs.size(); k++)
            {
                check_arr = (locs[j] - locs[i]) == loc_diffs[k];
                bool i_neigh = true;
                for(unsigned int l = 0; l < check_arr.size(); l++)
                {
                    i_neigh &= check_arr[l];
                }

                if(i_neigh)
                {
                    neighbours[i].push_back(j);
                    neigh_choice[i].push_back(k);
                    nn++;
                }
            }
        }
    }
}

void particle::field::field_type::gen_rand()
{
    if(ising)
    {
        testspin[0] = st_rand_int.gen()*2 - 1;
    }
    else
    {
        double phi = st_rand_double.gen()*2*M_PI;
        double cthet = 2*st_rand_double.gen()-1;
        double sthet = pow(1 - pow(cthet, 2), 0.5);
        testspin[0] = cos(phi)*sthet;
        testspin[1] = sin(phi)*sthet;
        testspin[2] = cthet;
    }
}

void particle::field::field_type::all_rand()
{
    for(unsigned int i = 0; i < spins.size(); i++)
    {
        this->gen_rand();
        this->set_rand(i);
    }
}

void particle::field::field_type::print_setup(const std::string filename,
    const std::string groupname,
    const int Tmax,
    const int Hmax)
{
    // Open existing file
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t f_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);

    // create datasets
    hsize_t *full_dims = (hsize_t*)malloc(sizeof(hsize_t)*(d+1));
    if (ising) {full_dims[d] = 1;}
    else {full_dims[d] = 3;}
    for(int i = 0; i < d; i++) {full_dims[i] = edgesize;}

    hid_t dspace_id = H5Screate_simple(d+1, full_dims, NULL);
    hid_t dset_id;
    std::stringstream nstream;
    std::string name;
    for(int i=0; i < Tmax; i++)
    {
        for(int j=0; j < Hmax; j++)
        {
            nstream << groupname << "/T_" << i << "-H_" << j;
            nstream >> name;
            nstream.clear();
            dset_id = H5Dcreate(f_id, name.c_str(), H5T_NATIVE_FLOAT,
                dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dclose(dset_id);
        }
    }

    free(full_dims);

    // close
    H5Fclose(f_id);
}

void particle::field::field_type::print(std::string filename,
    std::string arrname)
{
    // Find size of new array
    int t_size = 1;
    for (int i=0; i < d; i++) {t_size *= edgesize;}

    // Copy to float array
    float* new_x = alloc_1darr<float>(t_size);
    float* new_y = alloc_1darr<float>(t_size);
    float* new_z = alloc_1darr<float>(t_size);
    #pragma omp simd
    for(int i = 0; i < t_size; i++)
    {
        new_x[i] = 0;
        new_y[i] = 0;
        new_z[i] = 0;
    }
    #pragma omp simd
    for(int i = 0; i < spins.size(); i++)
    {
        int index = 0;
        for (int j = 0; j < d; j++)
        {
            index *= edgesize;
            index += locs[i][j];
        }
        new_x[index] = spins[i][0];
        if (!ising)
        {
            new_y[index] = spins[i][1];
            new_z[index] = spins[i][2];
        }
    }

    // Open existing file
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t f_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);

    // Open dataset
    hid_t dset_id = H5Dopen1(f_id, arrname.c_str());

    // Get slab and slice space
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    hsize_t *count = (hsize_t*)malloc(sizeof(hsize_t)*(d+1));
    hsize_t *offset = (hsize_t*)malloc(sizeof(hsize_t)*(d+1));
    count[d] = 1;
    offset[d] = 0;
    for(int i = 0; i < d; i++)
    {
        count[i] = edgesize;
        offset[i] = 0;
    }
    hid_t slice_space_id = H5Screate_simple(d+1, count, NULL);
    // for x
    hid_t slab_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(slab_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    H5Dwrite(dset_id, H5T_NATIVE_FLOAT, slice_space_id, slab_id, plist_id,
        new_x);
    H5Sclose(slab_id);

    if(!ising)
    {
        // for y
        offset[d] = 1;
        slab_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(slab_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Dwrite(dset_id, H5T_NATIVE_FLOAT, slice_space_id, slab_id, plist_id,
            new_y);
        H5Sclose(slab_id);
        // for z
        offset[d] = 2;
        slab_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(slab_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Dwrite(dset_id, H5T_NATIVE_FLOAT, slice_space_id, slab_id, plist_id,
            new_z);
        H5Sclose(slab_id);
    }

    // close
    H5Pclose(plist_id);
    H5Dclose(dset_id);
    H5Fclose(f_id);

    dealloc_1darr<float>(new_x);
    dealloc_1darr<float>(new_y);
    dealloc_1darr<float>(new_z);
    free(count);
    free(offset);
}

void particle::field::field_type::send_data(int dest_rank)
{
    // Send the metadata to set up recieves
    int metadata[6] =
        {spins.size(), d, ising, edgesize, periodic, loc_diffs.size()};
    MPI_Ssend(metadata, 6, MPI_INT, dest_rank, 0, MPI_COMM_WORLD);

    // Send the spins
    int nspinmod;
    if (ising) {nspinmod = 1;}
    else {nspinmod = 3;}
    double* spins_out = alloc_1darr<double>(spins.size()*nspinmod);
    int ind = 0;
    for(int i = 0; i < spins.size(); i++)
    {
        for(int j = 0; j < nspinmod; j++)
        {
            spins_out[ind] = spins[i][j];
            ind++;
        }
    }
    MPI_Ssend(spins_out, spins.size()*nspinmod, MPI_DOUBLE, dest_rank, 0,
        MPI_COMM_WORLD);

    // Send the locations
    int* loc_out = alloc_1darr<int>(spins.size()*d);
    ind = 0;
    for(int i = 0; i < spins.size(); i++)
    {
        for(int j = 0; j < d; j++)
        {
            loc_out[ind] = locs[i][j];
            ind++;
        }
    }
    MPI_Ssend(loc_out, spins.size()*d, MPI_INT, dest_rank, 0, MPI_COMM_WORLD);

    // Send the location diffs
    int* ld_out = alloc_1darr<int>(loc_diffs.size()*d);
    ind = 0;
    for(int i = 0; i < loc_diffs.size(); i++)
    {
        for(int j = 0; j < d; j++)
        {
            ld_out[ind] = loc_diffs[i][j];
            ind++;
        }
    }
    MPI_Ssend(ld_out, loc_diffs.size()*d, MPI_INT, dest_rank, 0,
        MPI_COMM_WORLD);

    // Send the Js
    double* J_out = alloc_1darr<double>(loc_diffs.size());
    for(int i = 0; i < loc_diffs.size(); i++)
    {
            J_out[i] = J_diffs[i];
    }
    MPI_Ssend(J_out, loc_diffs.size()*d, MPI_INT, dest_rank, 0,
        MPI_COMM_WORLD);

    dealloc_1darr<double>(spins_out);
    dealloc_1darr<int>(loc_out);
    dealloc_1darr<int>(ld_out);
    dealloc_1darr<double>(J_out);
}

void particle::field::field_type::recv_data(int src_rank)
{
    // Recv the metadata to set up other recieves
    int metadata[6];
    MPI_Status stat;
    MPI_Recv(metadata, 6, MPI_INT, src_rank, 0, MPI_COMM_WORLD, &stat);
    spins.resize(metadata[0]);
    locs.resize(metadata[0]);
    d = metadata[1];
    ising = metadata[2];
    edgesize = metadata[3];
    periodic = metadata[4];
    loc_diffs.resize(metadata[5]);
    J_diffs.resize(metadata[5]);

    this->set_default_spins();

    // Recv the spins
    int nspinmod;
    if (ising) {nspinmod = 1;}
    else {nspinmod = 3;}
    double* spins_out = alloc_1darr<double>(spins.size()*nspinmod);
    MPI_Recv(spins_out, spins.size()*nspinmod, MPI_DOUBLE, src_rank, 0,
        MPI_COMM_WORLD, &stat);
    int ind = 0;
    for(int i = 0; i < spins.size(); i++)
    {
        spins[i] = upspin;
        for(int j = 0; j < nspinmod; j++)
        {
            spins[i][j] = spins_out[ind];
            ind++;
        }
    }

    // Send the locations
    int* loc_out = alloc_1darr<int>(spins.size()*d);
    MPI_Recv(loc_out, spins.size()*d, MPI_INT, src_rank, 0, MPI_COMM_WORLD,
        &stat);
    ind = 0;
    for(int i = 0; i < spins.size(); i++)
    {
        locs[i] = blankloc;
        for(int j = 0; j < d; j++)
        {
            locs[i][j] = loc_out[ind];
            ind++;
        }
    }

    // Send the location diffs
    int* ld_out = alloc_1darr<int>(loc_diffs.size()*d);
    MPI_Recv(ld_out, loc_diffs.size()*d, MPI_INT, src_rank, 0,
        MPI_COMM_WORLD, &stat);
    ind = 0;
    for(int i = 0; i < loc_diffs.size(); i++)
    {
        locs[i] = blankloc;
        for(int j = 0; j < d; j++)
        {
            loc_diffs[i][j] = ld_out[ind];
            ind++;
        }
    }

    // Send the Js
    double* J_out = alloc_1darr<double>(loc_diffs.size());
    MPI_Recv(J_out, loc_diffs.size()*d, MPI_INT, src_rank, 0,
        MPI_COMM_WORLD, &stat);
    for(int i = 0; i < loc_diffs.size(); i++)
    {
            J_diffs[i] = J_out[i];
    }

    dealloc_1darr<double>(spins_out);
    dealloc_1darr<int>(loc_out);
    dealloc_1darr<int>(ld_out);
    dealloc_1darr<double>(J_out);

    // Recalculate the neighbours
    this->set_neigh();
}
