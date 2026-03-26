#pragma once
// Complex energy grid construction for SCF-MT calculation.
// Converted from src/POT/grids.f
//
// Makes a grid in the complex E-plane with three regions:
// 1. Near ecv: imaginary part increases quadratically
// 2. Horizontal: from ecv to xmu at constant imaginary part
// 3. Near xmu: imaginary part increases quadratically above real axis

#include <feff/types.hpp>

namespace feff::pot {

/// Construct complex energy grid for SCMT calculation.
///
/// @param ecv    Core-valence separation energy [Hartrees]
/// @param xmu    Fermi level [Hartrees]
/// @param negx   Maximum number of energy points (array dimension)
/// @param neg    Output: actual number of energy points
/// @param emg    Output: complex energy grid [negx]
/// @param step   Output: integration steps for Fermi level search [nflrx]
/// @param nflrx  Size of step array (controls grid density)
void grids(double ecv, double xmu, int negx, int& neg,
           FeffComplex* emg, double* step, int nflrx);

} // namespace feff::pot
