#pragma once
// Overlap densities from neighbors.
// Converted from src/POT/ovrlp.f
//
// Overlaps Coulomb potentials and electron densities for a given
// unique potential using Louck's spherical summation (sumax).

#include <feff/dimensions.hpp>

namespace feff::pot {

/// Overlap Coulomb potentials and electron densities for potential iph.
///
/// @param iph      Index of current unique potential (0-based)
/// @param iphat    Potential type per atom [natx] (0-based atom index)
/// @param rat      Atomic coordinates [3][natx] (column-major)
/// @param iatph    Representative atom for each potential type [nphx+1]
/// @param novr     Number of explicit overlap neighbors [nphx+1]
/// @param iphovr   Potential type of overlap neighbors [novrx][nphx+1] (flat)
/// @param nnovr    Number of each overlap neighbor [novrx][nphx+1] (flat)
/// @param rovr     Distance to overlap neighbors [novrx][nphx+1] (flat)
/// @param iz       Atomic number per potential [nphx+1]
/// @param nat      Number of atoms
/// @param rho      Free-atom density [251][nphx+2] (flat)
/// @param dmag     Spin density ratio [251][nphx+2] (flat, modified)
/// @param rhoval   Free-atom valence density [251][nphx+2] (flat)
/// @param vcoul    Free-atom Coulomb potential [251][nphx+2] (flat)
/// @param edens    Overlapped electron density [251][nphx+1] (flat, output)
/// @param edenvl   Overlapped valence density [251][nphx+1] (flat, output)
/// @param vclap    Overlapped Coulomb potential [251][nphx+1] (flat, output)
/// @param rnrm     Norman radii [nphx+1] (output for this iph)
void ovrlp(int iph, const int* iphat, const double* rat, const int* iatph,
           const int* novr, const int* iphovr, const int* nnovr, const double* rovr,
           const int* iz, int nat, const double* rho, double* dmag,
           const double* rhoval, const double* vcoul,
           double* edens, double* edenvl, double* vclap, double* rnrm);

} // namespace feff::pot
