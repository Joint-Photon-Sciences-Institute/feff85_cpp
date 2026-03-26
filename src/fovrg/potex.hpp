#pragma once
// Exchange potential for photoelectron in FOVRG.
// Converted from: src/FOVRG/potex.f

#include <feff/types.hpp>
#include <feff/dimensions.hpp>
#include "../atom/atom_types.hpp"

namespace feff::fovrg {

using feff::atom::FovrgState;

/// Calculate exchange potential for the photoelectron orbital.
/// Uses angular coefficients afgkc and radial integrals yzkrdc.
///
/// @param ps    Photoelectron large component (size nrptx)
/// @param qs    Photoelectron small component (size nrptx)
/// @param aps   Development coefficients for ps (size 10)
/// @param aqs   Development coefficients for qs (size 10)
/// @param jri   First interstitial grid point (0-based)
/// @param state FOVRG state (output: work.eg, work.ep, work.ceg, work.cep)
void potex(const FeffComplex ps[], const FeffComplex qs[],
           const FeffComplex aps[10], const FeffComplex aqs[10],
           int jri, FovrgState& state);

} // namespace feff::fovrg
