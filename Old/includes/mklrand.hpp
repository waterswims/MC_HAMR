#ifndef _MKLRAND
#define _MKLRAND

#include "mkl.h"

namespace mklrand
{
	///////////////////////////////////////////////////////////////////////////
	/// Base class for the random number generators.
	///////////////////////////////////////////////////////////////////////////
	class mkl_randbase
	{
	protected:
		int arr_size;
		int curr;
		VSLStreamStatePtr stream;
	public:
		////////////////////////////////////////////////////////////////////////
		/// Default constructor.
		////////////////////////////////////////////////////////////////////////
		mkl_randbase(){}
		////////////////////////////////////////////////////////////////////////
		/// Default destructor.
		////////////////////////////////////////////////////////////////////////
		~mkl_randbase(){vslDeleteStream(&stream);}
		////////////////////////////////////////////////////////////////////////
		/// Fill the buffer with new random numbers.
		////////////////////////////////////////////////////////////////////////
		virtual void fill(){}
		////////////////////////////////////////////////////////////////////////
		/// \brief Change the current random seed.
		///
		/// \param seed The new seed.
		////////////////////////////////////////////////////////////////////////
		void change_seed(int seed);
		////////////////////////////////////////////////////////////////////////
		/// \brief Save the state of the random number generator
		///
		/// \param name The name of the file which the state will be saved to.
		////////////////////////////////////////////////////////////////////////
		void save(const char* name) {vslSaveStreamF(stream, name);}
		////////////////////////////////////////////////////////////////////////
		/// \brief Load a random number generator state.
		///
		/// \param name The name of the file which the state loaded from.
		////////////////////////////////////////////////////////////////////////
		void load(const char* name) {vslDeleteStream(&stream); vslLoadStreamF(&stream, name);}
	};

	///////////////////////////////////////////////////////////////////////////
	/// Generator for uniform numbers between 0 and 1
	///////////////////////////////////////////////////////////////////////////
	class mkl_drand: public mkl_randbase
	{
	private:
		double *randarr;

	public:
		////////////////////////////////////////////////////////////////////////
		/// \brief Constructor
		///
		/// \param seed The inital seed of the random number generator.
		/// \param size The size of the buffer for RNG storage.
		////////////////////////////////////////////////////////////////////////
		mkl_drand(int size, int seed=1);
		////////////////////////////////////////////////////////////////////////
		/// Default destructor.
		////////////////////////////////////////////////////////////////////////
		~mkl_drand();
		////////////////////////////////////////////////////////////////////////
		/// Return a single random number.
		////////////////////////////////////////////////////////////////////////
		double gen();
		////////////////////////////////////////////////////////////////////////
		/// Fill the buffer with new random numbers.
		////////////////////////////////////////////////////////////////////////
		void fill();
	};

	///////////////////////////////////////////////////////////////////////////
	/// Generator for uniform integers between 0 and 1
	///////////////////////////////////////////////////////////////////////////
	class mkl_irand: public mkl_randbase
	{
	private:
		int *randarr;

	public:
		////////////////////////////////////////////////////////////////////////
		/// \brief Constructor
		///
		/// \param seed The inital seed of the random number generator.
		/// \param size The size of the buffer for RNG storage.
		////////////////////////////////////////////////////////////////////////
		mkl_irand(int size, int seed=1);
		////////////////////////////////////////////////////////////////////////
		/// Default destructor.
		////////////////////////////////////////////////////////////////////////
		~mkl_irand();
		////////////////////////////////////////////////////////////////////////
		/// Return a single random number.
		////////////////////////////////////////////////////////////////////////
		int gen();
		////////////////////////////////////////////////////////////////////////
		/// Fill the buffer with new random numbers.
		////////////////////////////////////////////////////////////////////////
		void fill();
	};

	///////////////////////////////////////////////////////////////////////////
	/// Generator for random numbers on a lognormal distribution
	///////////////////////////////////////////////////////////////////////////
	class mkl_lnrand: public mkl_randbase
	{
	private:
		double *randarr;
		double lmean, lsd;

	public:
		////////////////////////////////////////////////////////////////////////
		/// \brief Constructor
		///
		/// \param seed The inital seed of the random number generator.
		/// \param size The size of the buffer for RNG storage.
		/// \param m The logarithmic mean of the distribution.
		/// \param sd The logarithmic standard devaiation of the distribution.
		////////////////////////////////////////////////////////////////////////
		mkl_lnrand(double m, double sd, int size, int seed=1);
		////////////////////////////////////////////////////////////////////////
		/// Default destructor.
		////////////////////////////////////////////////////////////////////////
		~mkl_lnrand();
		////////////////////////////////////////////////////////////////////////
		/// Return a single random number.
		////////////////////////////////////////////////////////////////////////
		double gen();
		////////////////////////////////////////////////////////////////////////
		/// Fill the buffer with new random numbers.
		////////////////////////////////////////////////////////////////////////
		void fill();
	};
}

#endif