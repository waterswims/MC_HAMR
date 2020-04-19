#ifndef _STDRAND
#define _STDRAND

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace stdrand
{
    ///////////////////////////////////////////////////////////////////////////
    /// Base class for the random number generators.
    ///////////////////////////////////////////////////////////////////////////
    class std_randbase
    {
    protected:
        gsl_rng* generator;
    public:
    	////////////////////////////////////////////////////////////////////////
    	/// Default constructor.
    	////////////////////////////////////////////////////////////////////////
    	std_randbase() {generator = gsl_rng_alloc(gsl_rng_mt19937);}
    	////////////////////////////////////////////////////////////////////////
    	/// Default destructor.
    	////////////////////////////////////////////////////////////////////////
    	~std_randbase() {gsl_rng_free(generator);}
    	////////////////////////////////////////////////////////////////////////
    	/// \brief Change the current random seed.
    	///
    	/// \param seed The new seed.
    	////////////////////////////////////////////////////////////////////////
    	void change_seed(int seed) {gsl_rng_set(generator, seed);}
        ////////////////////////////////////////////////////////////////////////
    	/// Return a single random number.
    	////////////////////////////////////////////////////////////////////////
    	virtual double gen(){return 0;}
        ////////////////////////////////////////////////////////////////////////
    	/// Return a random vector.
    	////////////////////////////////////////////////////////////////////////
    	virtual void gen3(double* a, double* b, double* c){}
        ////////////////////////////////////////////////////////////////////////
    	/// \brief Jump forward 1e6 places
    	////////////////////////////////////////////////////////////////////////
    	void jump()
        {
            for(int i=0; i < 1e6; i++) {this->gen();}
        }
    };

    ///////////////////////////////////////////////////////////////////////////
    /// Generator for uniform numbers between 0 and 1
    ///////////////////////////////////////////////////////////////////////////
    class std_d_unirand : public std_randbase
    {
    public:
        ////////////////////////////////////////////////////////////////////////
    	/// \brief Constructor
    	///
    	/// \param seed The inital seed of the random number generator.
    	////////////////////////////////////////////////////////////////////////
        std_d_unirand(int seed);
        ////////////////////////////////////////////////////////////////////////
    	/// Default destructor.
    	////////////////////////////////////////////////////////////////////////
    	~std_d_unirand(){}
        ////////////////////////////////////////////////////////////////////////
    	/// Return a single random number.
    	////////////////////////////////////////////////////////////////////////
    	double gen(){return gsl_rng_uniform(generator);}
    };

    ///////////////////////////////////////////////////////////////////////////
    /// Generator for uniform integers between 0 and 1
    ///////////////////////////////////////////////////////////////////////////
    class std_i_unirand : public std_randbase
    {
    public:
        ////////////////////////////////////////////////////////////////////////
    	/// \brief Constructor
    	///
    	/// \param seed The inital seed of the random number generator.
    	////////////////////////////////////////////////////////////////////////
        std_i_unirand(int seed);
        ////////////////////////////////////////////////////////////////////////
    	/// Default destructor.
    	////////////////////////////////////////////////////////////////////////
    	~std_i_unirand(){}
        ////////////////////////////////////////////////////////////////////////
    	/// Return a single random number.
    	////////////////////////////////////////////////////////////////////////
    	double gen(){return gsl_rng_uniform_int(generator, 2);}
    };

    ///////////////////////////////////////////////////////////////////////////
    /// Generator for random numbers on a normal distribution
    ///////////////////////////////////////////////////////////////////////////
    class std_normrand : public std_randbase
    {
    private:
        double M, SD;
    public:
        ////////////////////////////////////////////////////////////////////////
    	/// \brief Constructor
    	///
    	/// \param m The mean of the normal distribution.
    	/// \param sdin The standard deviation of the normal distribution.
    	/// \param seed The inital seed of the random number generator. Defaults
    	///             to 1.
    	////////////////////////////////////////////////////////////////////////
        std_normrand(double m, double sdin, int seed=1);
        ////////////////////////////////////////////////////////////////////////
    	/// Default destructor.
    	////////////////////////////////////////////////////////////////////////
    	~std_normrand(){}
        ////////////////////////////////////////////////////////////////////////
    	/// Return a single random number.
    	////////////////////////////////////////////////////////////////////////
    	double gen(){return gsl_ran_gaussian(generator, SD) + M;}
    };

    ///////////////////////////////////////////////////////////////////////////
    /// Generator for random numbers on a lognormal distribution
    ///////////////////////////////////////////////////////////////////////////
    class std_lognormrand : public std_randbase
    {
    private:
        double M, SD;
    public:
        ////////////////////////////////////////////////////////////////////////
    	/// \brief Constructor
    	///
    	/// \param m The logarithmic mean of the lognormal distribution.
    	/// \param sdin The logarithmic standard deviation of the lognormal
        ///             distribution.
    	/// \param seed The inital seed of the random number generator.
    	////////////////////////////////////////////////////////////////////////
        std_lognormrand(double m, double sdin, int seed);
        ////////////////////////////////////////////////////////////////////////
    	/// Default destructor.
    	////////////////////////////////////////////////////////////////////////
    	~std_lognormrand(){}
        ////////////////////////////////////////////////////////////////////////
    	/// Return a single random number.
    	////////////////////////////////////////////////////////////////////////
    	double gen(){return gsl_ran_lognormal(generator, M, SD);}
    };

    ///////////////////////////////////////////////////////////////////////////
    /// Generator for random points on a sphere
    ///////////////////////////////////////////////////////////////////////////
    class std_sphere : public std_randbase
    {
    public:
        ////////////////////////////////////////////////////////////////////////
    	/// \brief Constructor
    	///
    	/// \param m The logarithmic mean of the lognormal distribution.
    	/// \param sdin The logarithmic standard deviation of the lognormal
        ///             distribution.
    	/// \param seed The inital seed of the random number generator.
    	////////////////////////////////////////////////////////////////////////
        std_sphere(int seed);
        ////////////////////////////////////////////////////////////////////////
    	/// Default destructor.
    	////////////////////////////////////////////////////////////////////////
    	~std_sphere(){}
        ////////////////////////////////////////////////////////////////////////
    	/// Return a random vector.
    	////////////////////////////////////////////////////////////////////////
    	void gen3(double* a, double* b, double* c)
            {gsl_ran_dir_3d(generator, a, b, c);}
    };
}

#endif