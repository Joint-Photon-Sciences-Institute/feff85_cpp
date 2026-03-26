#pragma once
// Make ri, beta (cos(beta)) path parameters for criterion calculations.
// Converted from: src/PATH/mrb.f

#include "path_data.hpp"

namespace feff::path {

/// Compute leg distances ri[0..npat] and cos(beta)[0..npat] for a path.
/// ipat[0..npat-1] are the atom indices in the path.
/// Uses shared AtomData for atom positions.
void mrb(int npat, const int ipat[], float ri[], float beta[],
         const AtomData& atoms);

/// Compute cos(angle) at point r between vectors (rm1->r) and (r->rp1).
/// Returns cos(beta). If degenerate (zero-length vector), returns 0.
float dotcos(const float rm1[3], const float r[3], const float rp1[3]);

} // namespace feff::path
