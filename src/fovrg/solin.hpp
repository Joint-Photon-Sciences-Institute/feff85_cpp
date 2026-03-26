#pragma once
// Inward solution of the Dirac equation (irregular solution).
// Converted from: src/FOVRG/solin.f
//
// Starts from spherical Hankel function outside muffin tin,
// uses flatv inward to iwkb, then Milne predictor-corrector to origin.

#include <feff/types.hpp>
#include <feff/dimensions.hpp>
#include "../atom/atom_types.hpp"

namespace feff::fovrg {

using feff::atom::DiracWorkspaceComplex;
using feff::atom::MeshParamsComplex;

/// Solve the Dirac equation inward from muffin-tin boundary.
///
/// @param en    Complex one-electron energy
/// @param fl    Power of first development term at origin
/// @param kap   Kappa quantum number
/// @param jri   First interstitial grid point (0-based)
/// @param imax  Last tabulation point (0-based)
/// @param ic3   Flag for c3 (spin-orbit) correction
/// @param vm    Spin-orbit correction potential (size nrptx)
/// @param iwkb  WKB switch point (0-based)
/// @param work  Dirac workspace (input: dv, eg, ep; output: gg, gp, ag, ap)
/// @param mesh  Mesh parameters
void solin(FeffComplex en, FeffComplex fl, int kap, int jri, int imax,
           int ic3, const FeffComplex vm[], int iwkb,
           DiracWorkspaceComplex& work, MeshParamsComplex& mesh);

} // namespace feff::fovrg
