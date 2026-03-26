#pragma once
// Spherical overlap summation (Louck's method).
// Converted from src/POT/sumax.f
//
// Performs eq 3.22 from Louck's "Augmented Plane Wave Method" (1967),
// using Simpson's rule on the logarithmic grid r(j) = exp(-8.8 + (j-1)*0.05).
// Linear interpolation is used at end caps.
//
// This version averages the contribution of one or more atoms of type 2
// at the location of atom 1.

namespace feff::pot {

/// Spherically summed overlap contribution (Louck's method).
///
/// @param rn     Distance from atom 1 to atom 2 (a.u.)
/// @param ann    Number of type 2 atoms (can be fractional)
/// @param aa2    Potential or density at atom 2 [250]
/// @param aasum  Output: spherically summed contribution added into
///               this array (can be called repeatedly to accumulate) [250]
void sumax(double rn, double ann, const double* aa2, double* aasum);

} // namespace feff::pot
