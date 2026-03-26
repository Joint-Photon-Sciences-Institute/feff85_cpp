#pragma once
// Heap path criterion: decides whether to add a path to the heap.
// Converted from: src/PATH/mcrith.f

#include "path_data.hpp"

namespace feff::path {

/// Compute heap importance factor xheap for a partial path.
/// Returns xheap = -1 if not defined (ss, triangles, or central atom ending).
/// npat: number of path atoms.
/// ipat[0..npat-1]: atom indices.
/// ri[0..npat]: leg distances.
/// indbet[0..npat]: nearest cos(beta) grid indices.
/// nncrit: number of criterion energy points.
/// fbetac: scattering amplitude array (flat, use fbetac_idx).
/// ckspc[0..nncrit-1]: |p| at criterion points.
void mcrith(int npat, const int ipat[], const float ri[], const int indbet[],
            const int ipot[], int nncrit,
            const float fbetac[], const float ckspc[],
            float& xheap);

} // namespace feff::path
