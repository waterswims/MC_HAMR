#include "../includes/protocol.hpp"

#include <cmath>
#include <hdf5.h>

void set_protocol(
    const int proto_code,
    float*& var1_list,
    float*& var2_list,
    int& var1_size,
    int& var2_size,
    int& var1_begin,
    int& var2_begin,
    int& var1_end,
    int& var2_end,
    int& var1_final,
    float* Hs,
    float* Ts,
    const int H_size,
    const int T_size)
{
    switch(proto_code)
    {
        case 1:
        var1_list = Hs;
        var2_list = Ts;
        var1_size = H_size;
        var2_size = T_size;
        var1_begin = var2_begin = 0;
        var1_end = H_size;
        var2_end = T_size;
        var1_final = H_size - 1;
        break;
        case 2:
        var1_list = Ts;
        var2_list = Hs;
        var1_size = T_size;
        var2_size = H_size;
        var1_begin = var2_begin = 0;
        var1_end = T_size;
        var2_end = H_size;
        var1_final = T_size - 1;
        break;
        case 4:
        var1_list = Ts;
        var2_list = Hs;
        var1_size = T_size;
        var2_size = H_size;
        var1_begin = 0;
        var2_end = -1;
        var1_end = T_size;
        var2_begin = H_size - 1;
        var1_final = T_size - 1;
        break;
        default:
        std::cout << "Invalid protocol, exiting..." << std::endl;
        exit(207);
    }
}

void incr_v1(
    const int proto_code,
    int& var1_curr)
{
    switch(proto_code)
    {
        case 1:
        case 2:
        case 4:
        var1_curr++;
        break;
        default:
        std::cout << "Invalid protocol, exiting..." << std::endl;
        exit(207);
    }
}

void incr_v2(
    const int proto_code,
    int& var2_curr)
{
    switch(proto_code)
    {
        case 1:
        case 2:
        var2_curr++;
        break;
        case 4:
        var2_curr--;
        break;
        default:
        std::cout << "Invalid protocol, exiting..." << std::endl;
        exit(207);
    }
}

bool check_rank_run(
    const int proto_code,
    const int i,
    const int comm_size,
    const int rank,
    const int var1_size)
{
    switch(proto_code)
    {
        case 1:
        case 2:
        case 4:
        return rank != i%comm_size;
        default:
        std::cout << "Invalid protocol, exiting..." << std::endl;
        exit(207);
    }
}

bool check_rank_latt(
    const int k,
    const int comm_size,
    const int rank,
    const int var1_size,
    int& sub_rank,
    int& sub_size,
    int& latt_rank,
    int& num_par)
{
    num_par = ceil(comm_size / float(var1_size));
    latt_rank = (rank / var1_size);
    sub_rank = rank-latt_rank*var1_size;
    sub_size = std::min<int>(var1_size, comm_size-latt_rank*var1_size);
    return latt_rank != k%num_par;
}

int count_sub_extra_opens(
    const int var1_size,
    const int var2_size,
    const int N_latts,
    const int rank,
    const int comm_size,
    const int k)
{
    int latt_prints = 2 * int(k == 0);

    int latt_rank = (rank / var1_size);
    int sub_rank = rank-latt_rank*var1_size;
    int num_in_sub = std::min<int>(var1_size, comm_size-latt_rank*var1_size);
    int num_v1s = (var1_size / num_in_sub)
        + int(sub_rank < var1_size%num_in_sub);

    int own_opens = (latt_prints + 1) * var2_size * num_v1s;

    int base_num_v1s = (var1_size / num_in_sub)
        + int(0 < var1_size%num_in_sub);

    int base_opens = (latt_prints + 1) * var2_size * base_num_v1s;

    return base_opens - own_opens;
}

int count_num_opens(
    const int var1_size,
    const int var2_size,
    const int N_latts,
    const int rank,
    const int comm_size)
{
    int latt_prints = 2 * int(rank == 0);

    int latt_rank = (rank / var1_size);
    int sub_rank = rank-latt_rank*var1_size;
    int num_in_sub = std::min<int>(var1_size, comm_size-latt_rank*var1_size);
    int num_v1s = (var1_size / num_in_sub) + int(sub_rank < var1_size%num_in_sub);

    int num_par = ceil(comm_size / float(var1_size));
    int num_reals = (N_latts / num_par) + int(latt_rank < N_latts%num_par);

    return (latt_prints + num_reals) * var2_size * num_v1s;
}

void extra_file_open(
    const std::string f_id)
{
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t file_id = H5Fopen(f_id.c_str(), H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);
    H5Fclose(file_id);
}
