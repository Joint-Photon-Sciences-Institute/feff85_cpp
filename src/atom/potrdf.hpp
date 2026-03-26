#pragma once
// Electron potential calculation for a given orbital.
// Converted from: src/ATOM/potrdf.f

#include "atom_types.hpp"

namespace feff::atom {

/// Calculate electron potential (Coulomb + exchange) for orbital ia.
/// Populates work.dv (direct potential), work.eg/ep (exchange), work.av/ceg/cep (dev coeffs).
/// Replaces Fortran potrdf(ia).
void potrdf(int ia, AtomState& state);

} // namespace feff::atom
