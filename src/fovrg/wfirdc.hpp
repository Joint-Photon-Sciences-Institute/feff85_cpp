#pragma once
// Calculate photoelectron orbital using LDA in the Dirac equation.
// Converted from: src/FOVRG/wfirdc.f

#include <feff/types.hpp>
#include <feff/dimensions.hpp>
#include "../atom/atom_types.hpp"

namespace feff::fovrg {

using feff::atom::FovrgState;

/// Calculate photoelectron orbital using LDA potential in Dirac equation.
///
/// @param eph   Complex one-electron energy of photoelectron
/// @param kap   Kappa quantum numbers for all orbitals (size 30, from config)
/// @param nmax  Last tabulation point for each orbital (size 30, from config)
/// @param vxc   Initial LDA potential for photoelectron (size nrptx)
/// @param ps    Output large component (size nrptx)
/// @param qs    Output small component (size nrptx)
/// @param aps   Output development coefficients for ps (size 10)
/// @param aqs   Output development coefficients for qs (size 10)
/// @param irr   Irregular solution flag (<0 = regular, >0 = irregular)
/// @param ic3   Flag for c3 (spin-orbit) correction
/// @param vm    Spin-orbit correction potential (size nrptx)
/// @param jri   First interstitial grid point (0-based)
/// @param iwkb  WKB switch point (0-based) [in/out]
/// @param state FOVRG state
void wfirdc(FeffComplex eph, int kap[], int nmax[], const FeffComplex vxc[],
            FeffComplex ps[], FeffComplex qs[], FeffComplex aps[10],
            FeffComplex aqs[10], int irr, int ic3, FeffComplex vm[],
            int jri, int& iwkb, FovrgState& state);

} // namespace feff::fovrg
