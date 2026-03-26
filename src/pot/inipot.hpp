#pragma once
// Initialize potential arrays to zero.
// Converted from src/POT/inipot.f

#include <feff/dimensions.hpp>

namespace feff::pot {

/// Zero out Dirac spinor components, valence density, valence potential,
/// and occupation numbers.
///
/// @param dgc     Large Dirac component [251][30][nphx+2] (flat, column-major)
/// @param dpc     Small Dirac component [251][30][nphx+2] (flat, column-major)
/// @param edenvl  Valence electron density [251][nphx+1]  (flat, column-major)
/// @param vvalgs  Valence ground-state potential [251][nphx+1] (flat)
/// @param xnmues  Occupation numbers [(lx+1)][nphx+1] (flat)
void inipot(double* dgc, double* dpc,
            double* edenvl, double* vvalgs, double* xnmues);

} // namespace feff::pot
