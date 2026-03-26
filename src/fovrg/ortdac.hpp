#pragma once
// Schmidt orthogonalization for FOVRG.
// Converted from: src/FOVRG/ortdac.f

#include <feff/types.hpp>
#include <feff/dimensions.hpp>
#include "../atom/atom_types.hpp"

namespace feff::fovrg {

using feff::atom::FovrgState;

/// Orthogonalize photoelectron orbital to all core orbitals of same symmetry.
/// Uses Schmidt orthogonalization procedure with dsordc overlap integrals.
///
/// @param ikap  Kappa quantum number for photoelectron
/// @param ps    Photoelectron large component (size nrptx) [in/out]
/// @param qs    Photoelectron small component (size nrptx) [in/out]
/// @param aps   Development coefficients for ps (size 10) [in/out]
/// @param aqs   Development coefficients for qs (size 10) [in/out]
/// @param state FOVRG state
void ortdac(int ikap, FeffComplex ps[], FeffComplex qs[],
            FeffComplex aps[10], FeffComplex aqs[10], FovrgState& state);

} // namespace feff::fovrg
