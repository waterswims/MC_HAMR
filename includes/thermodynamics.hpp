#ifndef _TDFUNCS
#define _TDFUNCS

#include "field_type.hpp"

#include <vector>
#include <functional>

namespace particle{ namespace td {
    class functionObject
    {
    private:
        std::function<double(field::field_type&,
            xt::xtensorf<double, xt::xshape<4>>&)> E_func;
        std::function<double(field::field_type&, int,
            xt::xtensorf<double, xt::xshape<4>>&)> dE_func;

    public:
        ////////////////////////////////////////////////////////////////////////
        /// Setup the thermodynamic functions
        ///
        /// \param useJ Should the exchange be included?
        /// \param useD Should the DMI be included?
        ////////////////////////////////////////////////////////////////////////
        void setup(bool useJ, bool useD);

        ////////////////////////////////////////////////////////////////////////
        /// Calculate the energy of a given lattice
        ///
        /// \param lattice The lattice
        /// \param H The external magnetic field
        /// \return The energy
        ////////////////////////////////////////////////////////////////////////
        double calc_E(field::field_type& lattice,
            xt::xtensorf<double, xt::xshape<4>>& H) {return E_func(lattice, H);}

        ////////////////////////////////////////////////////////////////////////
        /// Calculate the energy change after a spin flip of a single spin
        ///
        /// \param lattice The lattice
        /// \param position The chosen spin
        /// \param H The external magnetic field
        /// \return The change in energy
        ////////////////////////////////////////////////////////////////////////
        double calc_dE(field::field_type& lattice,
            int position, xt::xtensorf<double, xt::xshape<4>>& H)
            {return dE_func(lattice, position, H);}

        ////////////////////////////////////////////////////////////////////////
        /// Calculate the magnetisation of a lattice
        ///
        /// \param lattice The lattice
        /// \return The magnetisation vecotr of the lattice
        ////////////////////////////////////////////////////////////////////////
        xt::xtensorf<double, xt::xshape<4>> calc_M(field::field_type& lattice);

        ////////////////////////////////////////////////////////////////////////
        /// Calculate the sublattice magnetisation
        ///
        /// \param lattice The lattice
        /// \param subnumber Choice of sublattice
        /// \return The magnetisation vector of the sublattice
        ////////////////////////////////////////////////////////////////////////
        xt::xtensorf<double, xt::xshape<4>> calc_subM(
            field::field_type& lattice,
            int subnumber);

        ////////////////////////////////////////////////////////////////////////
        /// Calculate the 4-fold rotation sublattice magnetisation
        ///
        /// \param lattice The lattice
        /// \return The magnetisation vector of the sublattice
        ////////////////////////////////////////////////////////////////////////
        xt::xtensorf<double, xt::xshape<4>> calc_sub4M(
            field::field_type& lattice);

        ////////////////////////////////////////////////////////////////////////
        /// Calculate the topological charge of a given lattice
        ///
        /// \param lattice The lattice
        /// \return A vector containing the topological charge of each of the
        ///         slicesof the lattice
        ////////////////////////////////////////////////////////////////////////
        std::vector<double> calc_TC(field::field_type& lattice);
    };

    ///////////////////////////////////////////////////////////////////////////
    /// Calculate the solid angle between three vectors
    ///
    /// \param s1 First vector
    /// \param s2 Second vector
    /// \param s3 Third vector
    /// \return The solid angle
    ///////////////////////////////////////////////////////////////////////////
    double solid_angle(const xt::xtensorf<double, xt::xshape<4>>& s1,
                    const xt::xtensorf<double, xt::xshape<4>>& s2,
                    const xt::xtensorf<double, xt::xshape<4>>& s3);
}}

#endif
