#ifndef _FIELD
#define _FIELD

#include <vector>
#include <string>

#define XTENSOR_USE_XSIMD
#include <xtensor/xfixed.hpp>

namespace particle{ namespace field {
    ///////////////////////////////////////////////////////////////////////////
    /// Base class for fields
    ///////////////////////////////////////////////////////////////////////////
    class field_type
    {
    private:
        std::vector<xt::xtensorf<double, xt::xshape<4> > > spins;
        std::vector<xt::xtensorf<int, xt::xshape<4> > > locs;
        std::vector<std::vector<int> > neighbours;
        std::vector<std::vector<int> > adj;
        std::vector<std::vector<int> > neigh_choice;

        std::vector<xt::xtensorf<int, xt::xshape<4> > > loc_diffs;
        std::vector<double> J_diffs;
        std::vector<xt::xtensorf<double, xt::xshape<4> > > D_vecs;

        bool J_on;
        bool D_on;
        bool ising;
        bool periodic;
        int d;
        int edgesize;
        xt::xtensorf<double, xt::xshape<4>> upspin;
        xt::xtensorf<double, xt::xshape<4>> downspin;
        xt::xtensorf<double, xt::xshape<4>> testspin;
        xt::xtensorf<int, xt::xshape<4>> blankloc;

    public:
        ////////////////////////////////////////////////////////////////////////
        /// Default constructor
        ////////////////////////////////////////////////////////////////////////
        field_type() {}

        ////////////////////////////////////////////////////////////////////////
        /// Constructor
        ///
        /// \param ising_in Is this an ising field?
        /// \param periodic_in Is the field periodic?
        /// \param d_in The dimension of the field
        /// \param edgesize_in The edgelength of the lattice
        /// \param J_mod A modifier for the strength of the exchange
        /// \param D_mod A modifier for the strength of the DMI
        /// \param J_filename The filename where the exchanges and DMIs are
        ///                   stored
        ////////////////////////////////////////////////////////////////////////
        field_type(bool ising_in,
            bool periodic_in,
            int d_in,
            int edgesize_in,
            double J_mod,
            double D_mod,
            std::string J_filename);

        ////////////////////////////////////////////////////////////////////////
        /// Destructor
        ////////////////////////////////////////////////////////////////////////
        ~field_type() {}

        ////////////////////////////////////////////////////////////////////////
        /// Set default spins
        ////////////////////////////////////////////////////////////////////////
        void set_default_spins();

        ////////////////////////////////////////////////////////////////////////
        /// Access an individual spin
        ///
        /// \param index The index of the spin site
        /// \return A reference to the spin value as a valarray
        ////////////////////////////////////////////////////////////////////////
        xt::xtensorf<double, xt::xshape<4>>& access(int index) {return spins[index];}

        ////////////////////////////////////////////////////////////////////////
        /// Get the neighbours of a spin site
        ///
        /// \param index The index of the spin site
        /// \return A vector containing the indices of the neighbours of the
        ///         chosen spin
        ////////////////////////////////////////////////////////////////////////
        std::vector<int>& get_neigh(int index) {return neighbours[index];}

        ////////////////////////////////////////////////////////////////////////
        /// Get the neighbours of a spin site
        ///
        /// \param index The index of the spin site
        /// \return A vector containing the indices of the neighbours of the
        ///         chosen spin
        ////////////////////////////////////////////////////////////////////////
        std::vector<int>& get_adj(int index) {return adj[index];}

        ////////////////////////////////////////////////////////////////////////
        /// Get the exchanges with the neighbours of a spin site
        ///
        /// \param i The index of the spin site
        /// \param j The index of the chosen neighbour
        /// \return The exchange between the two spins
        ////////////////////////////////////////////////////////////////////////
        double& get_J(int i, int j) {return J_diffs[neigh_choice[i][j]];}

        ////////////////////////////////////////////////////////////////////////
        /// Get the DMI vector with the neighbours of a spin site
        ///
        /// \param i The index of the spin site
        /// \param j The index of the chosen neighbour
        /// \return The DMI vector between the two spins
        ////////////////////////////////////////////////////////////////////////
        xt::xtensorf<double, xt::xshape<4>>& get_D_vec(int i, int j)
            {return D_vecs[neigh_choice[i][j]];}

        ////////////////////////////////////////////////////////////////////////
        /// Get the location of a spin site
        ///
        /// \param index The index of the spin site
        /// \return A vector containing the location of the neighbours of the
        ///         chosen spin
        ////////////////////////////////////////////////////////////////////////
        xt::xtensorf<int, xt::xshape<4>>& get_loc(int index)
            {return locs[index];}

        ////////////////////////////////////////////////////////////////////////
        /// Add a new spin to the field of spins
        ///
        /// \param loc The location of the new spin
        ////////////////////////////////////////////////////////////////////////
        void add_spin(xt::xtensorf<int, xt::xshape<4>>& loc);

        ////////////////////////////////////////////////////////////////////////
        /// Determine the neighbours of the spins
        ////////////////////////////////////////////////////////////////////////
        void set_neigh();

        ////////////////////////////////////////////////////////////////////////
        /// Get the number of spins in the field
        ///
        /// \return The number of spins in the field
        ////////////////////////////////////////////////////////////////////////
        unsigned int get_size() {return spins.size();}

        ////////////////////////////////////////////////////////////////////////
        /// Get the dimension in the field
        ///
        /// \return The dimension of the field
        ////////////////////////////////////////////////////////////////////////
        int get_dim() {return d;}

        ////////////////////////////////////////////////////////////////////////
        /// Get the edgesize of the field
        ///
        /// \return The edgesize of the field
        ////////////////////////////////////////////////////////////////////////
        int get_edge() {return edgesize;}

        ////////////////////////////////////////////////////////////////////////
        /// Set a spin to the already generated random state
        ///
        /// \param index The location of the spin to be changed
        ////////////////////////////////////////////////////////////////////////
        void set_rand(int index) {spins[index] = testspin;}

        ////////////////////////////////////////////////////////////////////////
        /// Set a spin to the down state
        ///
        /// \param index The location of the spin to be changed
        ////////////////////////////////////////////////////////////////////////
        void set_up(int index) {spins[index] = upspin;}

        ////////////////////////////////////////////////////////////////////////
        /// Set a spin to the down state
        ///
        /// \param index The location of the spin to be changed
        ////////////////////////////////////////////////////////////////////////
        void set_down(int index) {spins[index] = downspin;}

        ////////////////////////////////////////////////////////////////////////
        /// Set a spin to the down state
        ///
        /// \param index The location of the spin to be changed
        ////////////////////////////////////////////////////////////////////////
        void set_spin(int index, xt::xtensorf<double, xt::xshape<4> >& in)
            {spins[index] = in;}

        ////////////////////////////////////////////////////////////////////////
        /// Generate a random spin state
        ////////////////////////////////////////////////////////////////////////
        void gen_rand();

        ////////////////////////////////////////////////////////////////////////
        /// Return the randomly generated state
        ///
        /// \return The random state as a valarray
        ////////////////////////////////////////////////////////////////////////
        xt::xtensorf<double, xt::xshape<4>>& get_rand() {return testspin;}

        ////////////////////////////////////////////////////////////////////////
        /// Set all spins to a random state
        ////////////////////////////////////////////////////////////////////////
        void all_rand();

        ////////////////////////////////////////////////////////////////////////
        /// Set all spins to a zero state
        ////////////////////////////////////////////////////////////////////////
        void all_zero();

        ////////////////////////////////////////////////////////////////////////
        /// Set up the HDF5 file for array printing
        ///
        /// /param filename The name of the HDF5 file
        /// /param groupname The name of the group within the file to store the
        ///                  lattices
        /// /param Tmax The number of temperatures
        /// /param Hmax The number of fields
        ////////////////////////////////////////////////////////////////////////
        void print_setup(const std::string filename,
            const std::string groupname,
            const int Tmax,
            const int Hmax);

        ////////////////////////////////////////////////////////////////////////
        /// Print the lattice to a file
        ///
        /// /param filename The name of the HDF5 file
        /// /param arrname The name of the particular array to be printed to
        ////////////////////////////////////////////////////////////////////////
        void print(std::string filename, std::string arrname);

        ////////////////////////////////////////////////////////////////////////
        /// Send the field to another process
        ///
        /// \param dest_rank The target process' rank
        ////////////////////////////////////////////////////////////////////////
        void send_data(int dest_rank);

        ////////////////////////////////////////////////////////////////////////
        /// Recieve the field from another process
        ///
        /// \param src_rank The rank of the process recieving from
        ////////////////////////////////////////////////////////////////////////
        void recv_data(int src_rank);

        ////////////////////////////////////////////////////////////////////////
        /// Are the exchanges on?
        ////////////////////////////////////////////////////////////////////////
        bool use_J() {return J_on;}

        ////////////////////////////////////////////////////////////////////////
        /// Are the dmis on?
        ////////////////////////////////////////////////////////////////////////
        bool use_D() {return D_on;}
    };

    ////////////////////////////////////////////////////////////////////////
    /// Find the short vector between two locations in periodic space
    ///
    /// \param A Position A
    /// \param B Position B
    /// \param C The output vector
    ////////////////////////////////////////////////////////////////////////
    void perioDiff(xt::xtensorf<int, xt::xshape<4>>& A,
        xt::xtensorf<int, xt::xshape<4>>& B,
        xt::xtensorf<int, xt::xshape<4>>& C);
}}

#endif
