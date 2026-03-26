#pragma once
// Non-diagonal Lagrange parameters for SCF calculation.
// Converted from: src/ATOM/lagdat.f
//
// Computes Lagrange multipliers eps(i,j) for pairs of orbitals with
// the same kappa quantum number but different occupation numbers.
// Used to enforce orthogonality between open and closed shells.

#include "atom_types.hpp"

namespace feff::atom {

/// Calculate non-diagonal Lagrange parameters.
///
/// @param ia   Orbital index (1-based). If positive, only Lagrange
///             parameters involving orbital ia are computed. If <= 0,
///             all Lagrange parameters are computed.
/// @param iex  Exchange flag. If 0, exchange contributions are omitted.
/// @param state Aggregate atom state
void lagdat(int ia, int iex, AtomState& state);

} // namespace feff::atom
