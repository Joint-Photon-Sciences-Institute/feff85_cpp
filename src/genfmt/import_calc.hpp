#pragma once

// Path importance factor calculation.
// Converted from GENFMT/import.f
//
// Named import_calc to avoid C++ keyword conflict with "import".

#include <feff/dimensions.hpp>
#include <feff/types.hpp>

namespace feff::genfmt {

/// Compute the importance factor of a path.
///
/// importance = deg * integral(|chi| * d|p|)
///
/// ne1: number of energy points on the main axis
/// nsp: number of spin channels
/// ik0: index of Fermi level on the energy grid (0-based)
/// reff: effective path length (half total)
/// deg: path degeneracy
/// ckmag: |ck| array
/// em: complex energy mesh
/// eref2: complex energy reference (spin-resolved)
/// cchi: complex chi array
/// xportx: reference importance (updated on first call if <= 0)
/// crit: output importance percentage (100 * xport / xportx)
void import_calc(int ne1, int nsp, int ik0, double reff, double deg,
                 const double ckmag[], const FeffComplex em[],
                 const FeffComplex eref2[][nspx], const FeffComplex cchi[],
                 double& xportx, double& crit);

} // namespace feff::genfmt
