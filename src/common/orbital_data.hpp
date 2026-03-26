#pragma once
// Orbital data and configuration for atomic elements.
// Converted from: src/COMMON/getorb.f
//
// Gets orbital data for chosen element: occupation numbers, quantum numbers,
// valence orbital information, and spin magnetization.

namespace feff::common {

/// Get orbital data for element with atomic number iz.
///
/// @param iz     Atomic number of desired element (1-99)
/// @param ihole  Index of core-hole orbital (0 = no hole)
/// @param xion   Ionicity (usually zero)
/// @param iunf   Unfolding flag for f-orbitals
/// @param norb   [out] Total number of orbitals
/// @param norbco [out] Number of core orbitals
/// @param iorb   [out] Index of orbital for making projections
/// @param iholep [out] Index of core hole orbital in compacted list
/// @param nqn    [out] Principal quantum number for each orbital (size 30)
/// @param nk     [out] Quantum number kappa for each orbital (size 30)
/// @param xnel   [out] Occupation for each orbital (size 30)
/// @param xnval  [out] Valence occupation for each orbital (size 30)
/// @param xmag   [out] Spin magnetization for each orbital (size 30)
void getorb(int iz, int ihole, double xion, int iunf,
            int& norb, int& norbco, int& iorb, int& iholep,
            int nqn[30], int nk[30], double xnel[30],
            double xnval[30], double xmag[30]);

} // namespace feff::common
