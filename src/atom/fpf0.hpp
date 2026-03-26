#pragma once
// Atomic scattering amplitude f0 and oscillator strengths.
// Converted from: src/ATOM/fpf0.f

#include "atom_types.hpp"
#include <feff/dimensions.hpp>

namespace feff::atom {

/// Calculate and output f0(Q) scattering amplitude and oscillator strengths.
/// Replaces Fortran fpf0().
void fpf0(int iz, int iholep, const double srho[251], const double dr[251],
           double hx, const double dgc0[251], const double dpc0[251],
           const double dgc[], const double dpc[],
           double eatom, const double xnel[30], int norb,
           const double eorb[30], const int kappa[30]);

} // namespace feff::atom
