#pragma once
// Diagnostic potential output.
// Converted from src/POT/wpot.f
//
// Writes potentials to potXX.dat files for each unique potential.
// Output includes free-atom and overlapped Coulomb potentials,
// densities, and total potentials.

#include <feff/dimensions.hpp>
#include <string>

namespace feff::pot {

/// Write potential data to potXX.dat files for diagnostics.
///
/// @param nph      Number of unique potentials
/// @param edens    Overlapped electron density [251][nphx+1] (flat)
/// @param imt      MT grid indices [nphx+1]
/// @param inrm     Norman grid indices [nphx+1]
/// @param rho      Free-atom density [251][nphx+2] (flat)
/// @param vclap    Overlapped Coulomb potential [251][nphx+1] (flat)
/// @param vcoul    Free-atom Coulomb potential [251][nphx+2] (flat)
/// @param vtot     Total potential [251][nphx+1] (flat)
/// @param ntitle   Number of title lines
/// @param title    Title lines [ntitle]
void wpot(int nph, const double* edens, const int* imt, const int* inrm,
          const double* rho, const double* vclap, const double* vcoul,
          const double* vtot, int ntitle, const std::string* title);

} // namespace feff::pot
