#pragma once
// Schmidt orthogonalization of orbitals.
// Converted from: src/ATOM/ortdat.f
//
// Orthogonalizes orbital ia to all orbitals of the same symmetry (kappa).
// If ia is positive, only orbital ia is orthogonalized.
// If ia is non-positive, all orbitals are orthogonalized sequentially.

#include "atom_types.hpp"

namespace feff::atom {

/// Schmidt orthogonalization.
///
/// @param ia    Orbital index (1-based). If positive, only orbital ia is
///              orthogonalized. If <= 0, all orbitals are orthogonalized.
/// @param state Aggregate atom state
void ortdat(int ia, AtomState& state);

} // namespace feff::atom
