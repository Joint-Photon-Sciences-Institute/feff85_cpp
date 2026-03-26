#pragma once
// Total energy calculation from orbital contributions.
// Converted from: src/ATOM/etotal.f

#include "atom_types.hpp"

namespace feff::atom {

/// Calculate total atomic energy from Coulomb, exchange, and Breit terms.
///   io: output file unit (0 = no output)
///   eatom: total energy in Hartrees (output)
/// Replaces Fortran etotal().
void etotal(int io, double& eatom, AtomState& state);

} // namespace feff::atom
