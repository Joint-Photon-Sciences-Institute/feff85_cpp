#pragma once
// Outward integration of the Dirac equation (complex version).
// Converted from: src/FOVRG/intout.f
//
// Uses Runge-Kutta start + Milne predictor-corrector.
// Solves the inhomogeneous Dirac equation:
//   p' - kap*p/r = -(en/cl - v)*g - eg/r
//   g' + kap*g/r = (2*cl + en/cl - v)*p + ep/r

#include <feff/types.hpp>
#include <feff/dimensions.hpp>
#include "../atom/atom_types.hpp"

namespace feff::fovrg {

using feff::atom::DiracWorkspaceComplex;
using feff::atom::MeshParamsComplex;

/// Integrate the inhomogeneous Dirac equation outward.
///
/// @param en    Complex one-electron energy (Hartrees)
/// @param i0    Starting grid point (0-based)
/// @param kap   Kappa quantum number
/// @param max0  Last point of tabulation (0-based)
/// @param ic3   Flag for c3 (spin-orbit) correction
/// @param vm    Spin-orbit correction potential (size nrptx)
/// @param work  Dirac workspace (gg, gp are input exchange / output wavefunctions)
/// @param mesh  Mesh parameters (hx, dr, ndor, np, idim)
void intout(FeffComplex en, int i0, int kap, int max0, int ic3,
            const FeffComplex vm[], DiracWorkspaceComplex& work,
            MeshParamsComplex& mesh);

} // namespace feff::fovrg
