#pragma once
// Edge energy lookup tables for all elements Z=1..98.
// Converted from src/XSPH/getedg.f
//
// Provides tabulated X-ray edge energies (in eV) from Elam tables
// for use in optical constant calculations.

#include <feff/dimensions.hpp>

namespace feff::xsph {

/// Look up edge energy from tabulated values.
///
/// If a valid entry exists for the given element and hole,
/// emu is set to the edge energy in Hartrees.
/// If the table entry is -1 (no data), emu is unchanged.
///
/// @param ihole  Core-hole index (1=K, 2=L1, 3=L2, 4=L3, ... up to 29)
/// @param iz     Atomic number (1..98)
/// @param emu    Edge energy in Hartrees (modified if valid entry found)
void getedg(int ihole, int iz, double& emu);

} // namespace feff::xsph
