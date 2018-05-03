#ifndef _PROT
#define _PROT

#include <iostream>

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
    const int T_size);

void incr_v1(
    const int proto_code,
    int& var1_curr);

void incr_v2(
    const int proto_code,
    int& var2_curr);

bool check_rank_run(
    const int proto_code,
    const int i,
    const int comm_size,
    const int rank,
    const int var1_size);

bool check_rank_latt(
    const int k,
    const int comm_size,
    const int rank,
    const int var1_size,
    int& sub_rank,
    int& sub_size,
    int& latt_rank,
    int& num_par);

int count_sub_extra_opens(
    const int var1_size,
    const int var2_size,
    const int N_latts,
    const int rank,
    const int comm_size,
    const int k);

int count_num_opens(
    const int var1_size,
    const int var2_size,
    const int N_latts,
    const int rank,
    const int comm_size);

void extra_file_open(
    const std::string f_id);

#endif
