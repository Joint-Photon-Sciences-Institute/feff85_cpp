#pragma once
// Angular Coulomb coefficients for the ATOM module.
// Converted from: muatco.f

#include "atom_types.hpp"

namespace feff::atom {

/// Compute angular coefficients for direct (Fk) and exchange (Gk) Coulomb
/// integrals using Wigner 3j symbols.
///
/// xnval[max_orb] — valence occupation numbers (0-based).
///                   If xnval[i] <= 0, orbital i is treated as core.
/// Results stored in ang.afgk[min(i,j)][max(i,j)][k/2] for Fk (direct)
///              and ang.afgk[max(i,j)][min(i,j)][k/2] for Gk (exchange).
void muatco(const double xnval[], ScfParams& scf,
            AngularCoefficients& ang, OrbitalConfig& config);

} // namespace feff::atom
