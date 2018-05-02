#include "../includes/array_alloc.hpp"
#include "../includes/state.hpp"
#include "../includes/functions.hpp"
#include "../includes/output.hpp"
#include "../includes/param_read.hpp"
#include "../includes/protocol.hpp"
#include "../includes/print_latt.hpp"

#ifdef __INTEL_COMPILER
#include "../includes/mklrand.hpp"
#define IRANDTYPE mklrand::mkl_irand
#define DRANDTYPE mklrand::mkl_drand
#define LNRANDTYPE mklrand::mkl_lnrand
#else
#include "../includes/stdrand.hpp"
#define IRANDTYPE stdrand::std_i_unirand
#define DRANDTYPE stdrand::std_d_unirand
#define LNRANDTYPE stdrand::std_lognormrand
#endif

#include <hdf5.h>
#include <cmath>

IRANDTYPE st_rand_int(1e7, 1);
DRANDTYPE st_rand_double(1e8, 1);

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
		std::cout << "Total number of processes is " << comm_size << std::endl;
	}

	// setup arguments
	int Samp_steps, Eq_steps, N_samp, protocol, N_latts;
	double J=1, K=0, k_b=1, beta=50, amean, asd, lmean, lsd, size;
	bool periodic, distributed, print_latt;
	char shape, hamil;
	std::string temp_name, field_name;
	read_all_vars(argv[1], size, J, k_b, periodic, shape, hamil, Samp_steps,
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
    LNRANDTYPE rand_ln(lmean, lsd, N_latts, 2*comm_size+rank);

    if(rank==0)
    {
		std::cout << num_Ts << " temperatures" << std::endl;
		std::cout << num_Hs << " field strengths" << std::endl;
        std::cout << N_samp << " samples per configuration" << std::endl;
    }

	// set protocol
	float* var1_list;
	float* var2_list;
	int v1_size, v2_size, v1_begin, v2_begin, v1_end, v2_end, v1_final;
	set_protocol(protocol, var1_list, var2_list, v1_size, v2_size, v1_begin,
		v2_begin, v1_end, v2_end, v1_final, Hs, Ts, num_Hs, num_Ts);

	// Start MC
	if (rank == 0)
	{
		std::cout << "Running MC..." << std::endl;
	}

	// Space for the average field
	field_type* summed_field = set_sum_latt(size, periodic, shape, hamil);
	int tc_size = summed_field->get_insize();

	// Create h5 file
	bool file_exists = false;
	bool** cpoint = alloc_2darr<bool>(v1_size, N_latts);
	std::string f_id = check_h5_file(J, size, amean, asd, k_b, K, shape,
									  hamil, periodic, protocol, distributed,
									  num_Ts, num_Hs, N_samp, N_latts, Ts, Hs,
									  v1_size, tc_size, cpoint, file_exists);
	if ((!file_exists) && print_latt)
	{
		summed_field->print_setup(f_id, "Av_Latt", num_Ts, num_Hs);
		summed_field->print_setup(f_id, "Sing_Latt", num_Ts, num_Hs);
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

	std::vector<double> mtemp;
	std::vector<double> tchargestemp;

	// if distrib, loop through latts
	for (int k = 0; k < N_latts; k++)
	{
		int sub_rank, sub_size, latt_rank, num_par;
		if (check_rank_latt(k, comm_size, rank, v1_size, sub_rank, sub_size,
							latt_rank, num_par))
		{
			int temp_int;
			double temp_double;
			for(int i = 0; i < 1e8; i++)
			{
				temp_int = st_rand_int.gen();
				temp_double = st_rand_double.gen();
				temp_double = rand_ln.gen();
			}
			continue;
		}
		if (distributed) {size = rand_ln.gen();}
		state base_state(size, periodic, shape, hamil, J, Hmax, k_b, Tmax, K,
			             args);
		state curr_state(base_state);
		nums = base_state.num_spins();
		s_nums = base_state.sub_num(0);

		base_state.equil(20*Eq_steps*nums);
		// main loop
		for (int i = v1_begin; i != v1_end && k < N_latts; incr_v1(protocol, i))
		{
			// Only carry on if to be run by this process
			if (check_rank_run(protocol, i, sub_size, sub_rank, v1_size))
			{
				continue;
			}

			// Grab the previous state from the previous rank unless first
			if (i != v1_begin)
			{
				int prev_rank = ((sub_rank+sub_size-1)%sub_size)
					+latt_rank*v1_size;
				base_state.recv_latt_data(prev_rank);
			}

			// Change Temp/field
			base_state.change_v1(protocol, var1_list[i]);
			// Equillibriation
			base_state.equil(Eq_steps*nums);

			// Send data to next rank, unless final rank
			if(sub_rank != sub_size-1 && i != v1_final)
			{
				int next_rank = ((sub_rank+1)%sub_size)+latt_rank*v1_size;
				base_state.send_latt_data(next_rank);
			}

			// Copy
			curr_state = base_state;

			// second loop, unless checkpointed
	        for (int j = v2_begin; j != v2_end && (!cpoint[i][k]);
				 incr_v2(protocol, j))
	        {
				// Change Temp/field
				curr_state.change_v2(protocol, var2_list[j]);
				// Equillibriation
				curr_state.equil(Eq_steps*nums);

				// Zero print field
				if(print_latt && k == 0)
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
					if(print_latt && k == 0)
					{
						curr_state.add_to_av(summed_field);
					}
				}

				// print values, lattice and checkpoint
				print_TD_h5(magx1, magy1, magz1, mag1, ener1, smagx1, smagy1,
					     smagz1, smag1, tcharges1, N_samp, protocol, i, j,
						 v2_size, tc_size, f_id, hamil, distributed, k);
				if (print_latt)
				{
					if (k == 0)
					{
						std::string lname = av_latt_name(protocol, i, j);
						summed_field->print(f_id, lname);
						lname = sing_latt_name(protocol, i, j);
						curr_state.ptf(f_id, lname);
					}
					else if (k/num_par == 0)
					{
						extra_file_open(f_id);
						extra_file_open(f_id);
					}
				}
	        }

			// Send data from final rank to first rank, unless final i
			if(sub_rank == sub_size-1 && i != v1_final)
			{
				int next_rank = (sub_rank+1)%sub_size+latt_rank*v1_size;
				base_state.send_latt_data(next_rank);
			}

		}
		// Extra HDF5 File open/closes for this k
		int ext_opens = count_sub_extra_opens(v1_size, v2_size, N_latts, rank,
			comm_size, k);
		for (int i = 0; i < ext_opens; i++)
		{
			extra_file_open(f_id);
		}
	}

	// Extra HDF5 File open/closes for whole sub set of ranks


	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0)
	{
		std::cout << std::endl;
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
