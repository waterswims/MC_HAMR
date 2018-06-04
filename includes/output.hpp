#ifndef _OUTP
#define _OUTP

#include "param_read.hpp"

#include <fstream>
#include <cstring>

void check_h5_file(stateOptions& stOpt,
                           simOptions& simOpt,
                           const int num_Ts,
                           const int num_Hs,
                           const float* Ts,
                           const float* Hs,
                           const int v1_size,
                           bool** checkp,
                           bool &file_exists);

void create_h5_file(stateOptions& stOpt,
                       simOptions& simOpt,
                       const int num_Ts,
                       const int num_Hs,
                       const float* Ts,
                       const float* Hs,
                       const int v1_size);

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
             stateOptions stOpt,
             simOptions simOpt,
             const int var1,
             const int var2,
             const int v2max,
             const int latt_num);

#endif
