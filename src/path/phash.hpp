#pragma once
// Path hashing for degeneracy checking.
// Converted from: src/PATH/phash.f

#include "path_data.hpp"

namespace feff::path {

/// Hash a path into a double-precision real number.
/// npat: number of path atoms. ipat[0..npat] (ipat[npat]=0 for central atom).
/// rx, ry, rz[0..npat-1]: path coordinates in standard frame.
/// dhash: output hash value.
void phash(int npat, const int ipat[], const float rx[], const float ry[],
           const float rz[], double& dhash, const AtomData& atoms);

} // namespace feff::path
