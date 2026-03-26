#pragma once
// Double-precision path parameter computation for output (ri, beta, eta).
// Converted from: src/PATH/mpprmd.f

#include "path_data.hpp"

namespace feff::path {

/// Compute double-precision path parameters for output to paths.dat.
/// npat: number of path atoms.
/// ipat[0..npat-1]: atom indices.
/// ri[0..npat]: leg distances (double precision).
/// beta[0..npat]: scattering angles in radians (double precision).
/// eta[0..npat]: Euler eta angles in radians (double precision).
void mpprmd(int npat, const int ipat[], double ri[], double beta[], double eta[],
            const AtomData& atoms);

} // namespace feff::path
