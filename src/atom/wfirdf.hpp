#pragma once
// Initial Thomas-Fermi orbital calculation via Dirac equation.
// Converted from: src/ATOM/wfirdf.f
//
// Calculates initial orbitals by:
// 1. Setting up the radial mesh with nucdev()
// 2. Computing Thomas-Fermi potential with dentfa()
// 3. Solving the Dirac equation for each orbital with soldir()

#include "atom_types.hpp"

namespace feff::atom {

/// Calculate initial orbitals from integration of the Dirac equation.
///
/// @param en    [in/out] One-electron energies (0-based, length max_orb)
/// @param ch    Ionicity (nuclear charge - number of electrons)
/// @param nq    Principal quantum numbers (0-based, length max_orb)
/// @param kap   Kappa quantum numbers (0-based, length max_orb)
/// @param nmax  [out] Number of tabulation points for each orbital (0-based)
/// @param ido   Option flag (only ido=1 is supported)
/// @param state Aggregate atom state
void wfirdf(double en[], double ch, int nq[], int kap[], int nmax[],
            int ido, AtomState& state);

} // namespace feff::atom
