#pragma once
// Coulomb potential development coefficients for FOVRG.
// Converted from: src/FOVRG/potdvp.f

#include <feff/types.hpp>
#include <feff/dimensions.hpp>
#include "../atom/atom_types.hpp"

namespace feff::fovrg {

using feff::atom::FovrgState;

/// Calculate potential development coefficients at origin.
/// Uses orbital density to compute Coulomb potential coefficients av.
///
/// @param state FOVRG state (reads orb, config, nuclear; writes work.av)
void potdvp(FovrgState& state);

} // namespace feff::fovrg
