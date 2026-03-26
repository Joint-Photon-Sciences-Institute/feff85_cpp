#pragma once

// Spherical wave factors c_l^(m) z^m / m!
// Converted from GENFMT/sclmz.f

#include "genfmt_data.hpp"

namespace feff::genfmt {

/// Compute CLM(Z) for current leg.
/// Makes clm(z) (eq B11 of Rehr-Albers).
///
/// clmi(il, im, ileg) is filled for ileg, elements clm(0,0) -> clm(lmax+1, mmax+1).
///
/// rho: complex rho array (rho[ileg] = ck * ri[ileg])
/// lmaxp1: lmax + 1 for this leg
/// mmaxp1: mmax + 1
/// ileg: 0-based leg index
void sclmz(const FeffComplex rho[], int lmaxp1, int mmaxp1, int ileg,
           ClmzData& clmz);

} // namespace feff::genfmt
