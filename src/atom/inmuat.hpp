#pragma once
// Initialize orbital occupations and SCF parameters.
// Converted from: src/ATOM/inmuat.f
//
// Sets up orbital data using getorb() from COMMON module, initializes
// Lagrange parameters, precision tests, and shell occupation flags.

#include "atom_types.hpp"

namespace feff::atom {

/// Initialize orbital occupations and SCF configuration.
///
/// @param ihole   Index of core-hole orbital (0 = no hole)
/// @param xionin  Ionicity
/// @param iunf    Unfolding flag for f-orbitals
/// @param xnval   [out] Valence occupation for each orbital (0-based, length max_orb)
/// @param iholep  [out] Index of core hole orbital in compacted list
/// @param xmag    [out] Spin magnetization for each orbital (0-based, length max_orb)
/// @param iorb    [out] Index of orbital for making projections
/// @param state   Aggregate atom state (config, scf, lagrange, nuclear, mesh populated)
void inmuat(int ihole, double xionin, int iunf, double xnval[],
            int& iholep, double xmag[], int& iorb, AtomState& state);

} // namespace feff::atom
