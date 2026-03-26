#pragma once
// Angular exchange coefficients for FOVRG.
// Converted from: src/FOVRG/muatcc.f
//
// Computes b_k(ikap,j) coefficients using Wigner 3j symbols.

#include <feff/types.hpp>
#include <feff/dimensions.hpp>
#include "../atom/atom_types.hpp"

namespace feff::fovrg {

using feff::atom::FovrgState;

/// Calculate angular exchange coefficients afgkc.
///
/// @param xnval Valence occupation numbers (size 30)
/// @param state FOVRG state (modifies angular coefficients)
void muatcc(const double xnval[30], FovrgState& state);

} // namespace feff::fovrg
