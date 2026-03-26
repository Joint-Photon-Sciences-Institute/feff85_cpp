#pragma once
// Outward solution of the Dirac equation with power-series start.
// Converted from: src/FOVRG/solout.f
//
// Computes development coefficients at origin using Desclaux expansion
// or c3-corrected expansion, then calls intout for numerical integration,
// then uses flatv for the region beyond iwkb.

#include <feff/types.hpp>
#include <feff/dimensions.hpp>
#include "../atom/atom_types.hpp"

namespace feff::fovrg {

using feff::atom::DiracWorkspaceComplex;
using feff::atom::MeshParamsComplex;

/// Solve the Dirac equation outward from origin.
///
/// @param en    Complex one-electron energy
/// @param fl    Power of first development term at origin
/// @param agi   Initial large component dev. coefficient
/// @param api   Initial small component dev. coefficient
/// @param kap   Kappa quantum number
/// @param jri   First interstitial grid point (0-based)
/// @param max0  Last tabulation point (0-based)
/// @param ic3   Flag for c3 (spin-orbit) correction
/// @param vm    Spin-orbit correction potential (size nrptx)
/// @param iwkb  WKB switch point (0-based)
/// @param work  Dirac workspace (input: dv, eg, ep; output: gg, gp, ag, ap)
/// @param mesh  Mesh parameters
void solout(FeffComplex en, FeffComplex fl, FeffComplex agi, FeffComplex api,
            int kap, int jri, int max0, int ic3, const FeffComplex vm[],
            int iwkb, DiracWorkspaceComplex& work, MeshParamsComplex& mesh);

} // namespace feff::fovrg
