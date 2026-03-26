#pragma once

// Many-pole self-energy grid calculation
// Converted from src/EXCH/mpse.f
//
// Calculates the many-pole self-energy at each energy point
// and a few r points, saving the results to mpse.bin for later
// use by xcpot.

#include <feff/types.hpp>
#include <feff/dimensions.hpp>

namespace feff::exch {

/// Many-pole self-energy grid calculation.
///
/// Reads exc.dat for pole parameters, calculates Sigma on a grid of
/// NRPts density-parameter points from RsMin to RsMax, and writes
/// the results to mpse.bin.
///
/// @param edens   Electron density array, edens[nphx+1][252] (row = potential, col = radial)
///                Fortran layout: edens(251, 0:nphx) -- 1-based radial, 0-based potential
/// @param jintrs  Index of last point before interstitial level (0-based)
/// @param e       Current energy (complex, a.u.)
/// @param ipl     Many-pole self energy control
/// @param mu      Fermi energy (a.u.)
void mpse(const double edens[][252], int jintrs, FeffComplex e,
          int ipl, double mu);

} // namespace feff::exch
