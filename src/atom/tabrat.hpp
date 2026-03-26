#pragma once
// Orbital statistics tabulation.
// Converted from: src/ATOM/tabrat.f

#include "atom_types.hpp"

namespace feff::atom {

/// Tabulate orbital results and statistics (average r^n values, overlaps).
/// Replaces Fortran tabrat().
void tabrat(AtomState& state);

} // namespace feff::atom
