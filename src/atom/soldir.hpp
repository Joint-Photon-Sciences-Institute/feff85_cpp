#pragma once
// Dirac equation solver with iterative energy adjustment.
// Converted from: src/ATOM/soldir.f
//
// Solves the Dirac equation for bound states by integrating with intdir()
// and iteratively adjusting the energy to match boundary conditions and
// achieve the correct number of nodes.

#include "atom_types.hpp"

namespace feff::atom {

/// Solve the Dirac equation for a bound-state orbital.
///
/// @param en      [in/out] One-electron energy (Hartrees, negative)
/// @param fl      Power of the first development term at the origin
/// @param agi     [in/out] First development coefficient (large component)
/// @param api     [in/out] First development coefficient (small component)
/// @param ainf    [in/out] Initial value for large component at boundary
/// @param nq      Principal quantum number
/// @param kap     Kappa quantum number
/// @param max0    [out] Last tabulation point index (1-based)
/// @param ifail   [out] Set to 1 if iteration limits exceeded (non-fatal warning)
/// @param work    Dirac workspace (dv, av, eg, ep, gg, gp, ag, ap arrays)
/// @param mesh    Radial mesh parameters
/// @param solver  Dirac solver state
/// @param error   Error state
void soldir(double& en, double fl, double& agi, double& api, double& ainf,
            int nq, int kap, int& max0, int& ifail,
            DiracWorkspaceReal& work, MeshParamsReal& mesh,
            DiracSolverState& solver, ErrorState& error);

} // namespace feff::atom
