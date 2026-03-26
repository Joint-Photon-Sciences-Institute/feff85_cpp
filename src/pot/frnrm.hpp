#pragma once
// Norman radius from electron density.
// Converted from src/POT/frnrm.f
//
// Finds the Norman sphere radius by integrating 4*pi*rho*r^2 outward
// until the integral equals Z (atomic number).  Uses extended Simpson
// integration with Newton-Raphson correction.

#include <feff/dimensions.hpp>

namespace feff::pot {

/// Find Norman radius from overlapped electron density.
///
/// @param rho   Overlapped density (4*pi*density) on Loucks grid [nrptx]
/// @param iz    Atomic number of the atom
/// @param rnrm  Output: Norman radius (a.u.)
void frnrm(const double* rho, int iz, double& rnrm);

} // namespace feff::pot
