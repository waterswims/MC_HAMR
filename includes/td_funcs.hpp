#ifndef _TDFUNCS
#define _TDFUNCS

#include "field_type.hpp"

#include <vector>
#include <valarray>

namespace particle{ namespace funcs {
    ///////////////////////////////////////////////////////////////////////////
    /// Calculate the energy of a given lattice
    ///
    /// \param lattice The lattice
    /// \param H The external magnetic field
    /// \return The energy
    ///////////////////////////////////////////////////////////////////////////
    double calc_E(field::field_type& lattice, std::valarray<double> H);

    ///////////////////////////////////////////////////////////////////////////
    /// Calculate the energy change after a spin flip of a single spin
    ///
    /// \param lattice The lattice
    /// \param position The chosen spin
    /// \param H The external magnetic field
    /// \return The change in energy
    ///////////////////////////////////////////////////////////////////////////
    double calc_dE(field::field_type& lattice, int position, std::valarray<double> H);

    ///////////////////////////////////////////////////////////////////////////
    /// Calculate the magnetisation of a lattice
    ///
    /// \param lattice The lattice
    /// \return The magnetisation vecotr of the lattice
    ///////////////////////////////////////////////////////////////////////////
    std::valarray<double> calc_M(field::field_type& lattice);

    ///////////////////////////////////////////////////////////////////////////
    /// Calculate the sublattice magnetisation
    ///
    /// \param lattice The lattice
    /// \param subnumber Choice of sublattice
    /// \return The magnetisation vecotr of the sublattice
    ///////////////////////////////////////////////////////////////////////////
    std::valarray<double> calc_subM(field::field_type& lattice, int subnumber);

    ///////////////////////////////////////////////////////////////////////////
    /// Calculate the topological charge of a given lattice
    ///
    /// \param lattice The lattice
    /// \return A vector containing the topological charge of each of the slices
    ///         of the lattice
    ///////////////////////////////////////////////////////////////////////////
    std::vector<double> calc_TC(field::field_type& lattice);

    ///////////////////////////////////////////////////////////////////////////
    /// Calculate the solid angle between three vectors
    ///
    /// \param s1 First vector
    /// \param s2 Second vector
    /// \param s3 Third vector
    /// \return The solid angle
    ///////////////////////////////////////////////////////////////////////////
    double solid_angle(const std::valarray<double> &s1,
                    const std::valarray<double> &s2,
                    const std::valarray<double> &s3);
}}

#endif
