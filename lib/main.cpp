#include "../includes/array_alloc.hpp"
#include "../includes/state.hpp"
#include "../includes/mklrand.hpp"
#include "../includes/functions.hpp"
#include "../includes/output.hpp"
#include "../includes/param_read.hpp"
#include "../includes/protocol.hpp"
#include "../includes/print_latt.hpp"

#include <iostream>
#include <vector>
#include <fstream>
#include <cstring>
#include <hdf5.h>

mkl_irand st_rand_int(1e7, 1);
mkl_drand st_rand_double(1e8, 1);

int main(int argc, char **argv)
{
    //Start MPI off
	MPI_Init(&argc, &argv);
	int rank, comm_size;
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;
  	if(rank==0)
	{
		cout << "Total number of processes is " << comm_size << endl;
	}

	// setup arguments
	int Samp_steps, Eq_steps, N_samp, protocol, N_latts;
	double J=1, K=0, k=1, beta=50, amean, asd, lmean, lsd, size;
	bool periodic, distributed, print_latt;
	char shape, hamil;
	string temp_name, field_name;
	read_all_vars(argv[1], size, J, k, periodic, shape, hamil, Samp_steps,
		N_samp, Eq_steps, N_latts, beta, distributed, amean, asd, temp_name,
		field_name, protocol, K, print_latt);
	double args[] = {beta};
	AtoLn(amean, asd, lmean, lsd);

	// load temperatures and fields
	float* Ts;
	int num_Ts;
	float* Hs;
	int num_Hs;
	load_Hs_Ts(temp_name, Ts, num_Ts, field_name, Hs, num_Hs);
	double Tmin(Ts[num_Ts-1]), Tmax(Ts[0]);
	double Hmin(Hs[num_Hs-1]), Hmax(Hs[0]);

	//Set the random seeds of the generators
	st_rand_int.change_seed(rank);
	st_rand_double.change_seed(comm_size+rank);

	// New Generator for Lognormal
    mkl_lnrand rand_ln(lmean, lsd, N_latts, 2*comm_size+rank);

    if(rank==0)
    {
		cout << num_Ts << " temperatures" << endl;
		cout << num_Hs << " field strengths" << endl;
        cout << N_samp << " samples per configuration" << endl;
    }

	// set protocol
	float* var1_list;
	float* var2_list;
	int v1_size, v2_size;
	set_protocol(protocol, var1_list, var2_list, v1_size, v2_size,Hs, Ts,
		num_Hs, num_Ts);

	// Start MC
	if (rank == 0)
	{
		cout << "Running MC..." << endl;
	}

	// Space for the average field
	field_type* summed_field = set_sum_latt(size, periodic, shape, hamil);
	int tc_size = summed_field->get_insize();

	// Create h5 file
	bool file_exists = false;
	bool** cpoint = alloc_2darr<bool>(v1_size, N_latts);
	std::string f_id = check_h5_file(J, size, amean, asd, k, K, shape,
									  hamil, periodic, protocol, distributed,
									  num_Ts, num_Hs, N_samp, N_latts, Ts, Hs,
									  v1_size, tc_size, cpoint, file_exists);
	if ((!file_exists))
	{
		summed_field->print_setup(f_id, num_Ts, num_Hs);
	}

	// Storage Variables
	int nums, s_nums;
	float* mag1 = alloc_1darr<float>(N_samp);
	float* ener1 = alloc_1darr<float>(N_samp);
	float* magx1 = alloc_1darr<float>(N_samp);
	float* magy1 = alloc_1darr<float>(N_samp);
	float* magz1 = alloc_1darr<float>(N_samp);
	float* smag1 = alloc_1darr<float>(N_samp);
	float* smagx1 = alloc_1darr<float>(N_samp);
	float* smagy1 = alloc_1darr<float>(N_samp);
	float* smagz1 = alloc_1darr<float>(N_samp);
	float** tcharges1 = alloc_2darr<float>(tc_size, N_samp);

	vector<double> mtemp;
	vector<double> tchargestemp;

	// if distrib, loop through latts
	for (int k = 0; k < N_latts; k++)
	{
		if (distributed) {size = rand_ln.gen();}
		state base_state(size, periodic, shape, hamil, J, Hmax, k, Tmax, K,
			             args);
		state curr_state(base_state);
		nums = base_state.num_spins();
		s_nums = base_state.sub_num(0);
		// main loop
		for (int i = 0; i < v1_size; i++)
		{
			// Change Temp/field
			base_state.change_v1(protocol, var1_list[i]);
			// Equillibriation
			base_state.equil(Eq_steps*nums);
			// Only carry on if to be run by this process and not checkpointed
			if (i%comm_size != rank || cpoint[i][k])
			{
				continue;
			}
			// Copy
			curr_state = base_state;
	        for (int j = 0; j < v2_size; j++)
	        {
				// Change Temp/field
				curr_state.change_v2(protocol, var2_list[j]);
				// Equillibriation
				curr_state.equil(Eq_steps*nums);

				// Zero print field
				if(print_latt)
				{
					summed_field->allzero();
				}

				// Get Samples
				for (int ns = 0; ns < N_samp; ns++)
				{
					curr_state.equil(Samp_steps*nums);
					mtemp = curr_state.magnetisation();
					if (hamil != 'i' && hamil != 'I')
					{
						magx1[ns] = mtemp[0];
						magy1[ns] = mtemp[1];
						magz1[ns] = mtemp[2];
					}
					else
					{
						magz1[ns] = mtemp[0];
					}

					mag1[ns] = norm(mtemp);
					ener1[ns] = curr_state.energy();

					mtemp = curr_state.submag(0);
					if (hamil != 'i' && hamil != 'I')
					{
						smagx1[ns] = mtemp[0];
						smagy1[ns] = mtemp[1];
						smagz1[ns] = mtemp[2];
					}
					else
					{
						smagz1[ns] = mtemp[0];
					}
					smag1[ns] = norm(mtemp);

					if (hamil != 'i' && hamil != 'I')
					{
						tchargestemp = curr_state.tcharge();
						for(int tc_idx=0; tc_idx < tc_size; tc_idx++)
						{
							tcharges1[tc_idx][ns] = tchargestemp[tc_idx];
						}
					}

					// Add sample to print lattice
					if(print_latt && !distributed)
					{
						curr_state.add_to_av(summed_field);
					}
				}
				// print values, lattice and checkpoint
				print_TD_h5(magx1, magy1, magz1, mag1, ener1, smagx1, smagy1,
					     smagz1, smag1, tcharges1, N_samp, protocol, i, j,
						 v2_size, tc_size, f_id, hamil, distributed, k);
				if (print_latt && !distributed)
				{
					std::string lname = latt_name(protocol, i, j);
					summed_field->print(f_id, lname);
				}
	        }
		}
		// Extra File open/closes
		if(rank >= v1_size - (v1_size/comm_size)*comm_size &&
			!cpoint[v1_size-1][k])
		{
			for(int i = 0; i < v2_size; i++)
			{
				hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
			    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
			    hid_t file_id = H5Fopen(f_id.c_str(), H5F_ACC_RDWR, plist_id);
			    H5Pclose(plist_id);
				H5Fclose(file_id);

				if (print_latt && !distributed)
				{
					hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
				    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
				    hid_t file_id = H5Fopen(f_id.c_str(), H5F_ACC_RDWR,
						plist_id);
				    H5Pclose(plist_id);
					H5Fclose(file_id);
				}
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0)
	{
		cout << endl;
	}

	// Destroy arrays
	dealloc_1darr<float>(Ts);
	dealloc_1darr<float>(Hs);
	dealloc_1darr<float>(mag1);
	dealloc_1darr<float>(ener1);
	dealloc_1darr<float>(magx1);
	dealloc_1darr<float>(magy1);
	dealloc_1darr<float>(magz1);
	dealloc_1darr<float>(smag1);
	dealloc_1darr<float>(smagx1);
	dealloc_1darr<float>(smagy1);
	dealloc_1darr<float>(smagz1);
	dealloc_2darr<bool>(v1_size, cpoint);
	dealloc_2darr<float>(tc_size, tcharges1);

    // Finish program
	MPI_Finalize();

    return 0;
}