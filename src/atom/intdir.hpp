#pragma once
// Dirac equation integrator using 5-point predictor-corrector method.
// Converted from: src/ATOM/intdir.f
//
// Solves the inhomogeneous Dirac equation by outward integration from
// the origin, inward integration from the boundary, and matching at a
// turning point.

#include "atom_types.hpp"

namespace feff::atom {

/// Integrate the (inhomogeneous) Dirac equation.
///
/// @param gg      [in/out] Large component: exchange terms on entry, wave function on exit (0-based, length atom_grid)
/// @param gp      [in/out] Small component: exchange terms on entry, wave function on exit
/// @param ag      [out] Development coefficients for gg at origin (length max_dev)
/// @param ap      [out] Development coefficients for gp at origin
/// @param ggmat   [out] Large component value at matching point (inward integration)
/// @param gpmat   [out] Small component value at matching point (inward integration)
/// @param en      One-electron energy (Hartrees, negative)
/// @param fl      Power of the first development term at the origin
/// @param agi     Initial value of the first development coefficient (large component)
/// @param api     Initial value of the first development coefficient (small component)
/// @param ainf    Initial value for large component at dr[max0-1]
/// @param max0    [in/out] Last tabulation point index (1-based Fortran convention preserved internally)
/// @param work    Dirac workspace (dv, av arrays = direct potential)
/// @param mesh    Radial mesh parameters
/// @param solver  Dirac solver state (ell, fk, ccl, imm, nd, node, mat)
/// @param error   Error state (numerr set on failure)
void intdir(double gg[], double gp[], double ag[], double ap[],
            double& ggmat, double& gpmat,
            double en, double fl, double agi, double api, double ainf,
            int& max0,
            DiracWorkspaceReal& work, MeshParamsReal& mesh,
            DiracSolverState& solver, ErrorState& error);

} // namespace feff::atom
