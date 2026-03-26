#pragma once
// Coulomb potential with charge normalization.
// Converted from src/POT/coulom.f
//
// Updates the Coulomb potential from the difference between new and old
// valence densities. Supports two normalization modes:
//   icoul=1: explicit charge normalization using neighbors
//   icoul=0: Norman picture normalization (default for infinite solids)

#include <feff/dimensions.hpp>

namespace feff::pot {

/// Update Coulomb potential from density change.
///
/// @param icoul    Normalization mode (0=Norman picture, 1=explicit)
/// @param npot     Number of potentials
/// @param ilast    Last grid index per potential [nphx+1]
/// @param rhoval   New valence density [251][nphx+2] (flat)
/// @param edenvl   Old valence density [251][nphx+1] (flat)
/// @param edens    Total electron density [251][nphx+1] (flat)
/// @param nat      Number of atoms
/// @param rat      Atomic coordinates [3][natx] (column-major)
/// @param iatph    Representative atom per potential [nphx+1]
/// @param iphat    Potential type per atom [natx]
/// @param rnrm     Norman radii [nphx+1]
/// @param dq       Charge transfer per potential [nphx+1]
/// @param iz       Atomic numbers [nphx+1]
/// @param vclap    Coulomb potential [251][nphx+1] (flat, modified)
void coulom(int icoul, int npot, const int* ilast, const double* rhoval,
            const double* edenvl, const double* edens,
            int nat, const double* rat, const int* iatph, const int* iphat,
            double* rnrm, const double* dq, const int* iz, double* vclap);

/// Helper function for Norman normalization: integral of
/// (4*pi*rho) * r^2 * (1/r0 - 1/r) where rho = a*r + b
double fab(double aa, double bb, double r0, double r1, double r2);

} // namespace feff::pot
