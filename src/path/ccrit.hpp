#pragma once
// Combined path criterion: heap + keep decisions.
// Converted from: src/PATH/ccrit.f

#include "path_data.hpp"

namespace feff::path {

/// Evaluate combined heap/keep criteria for a path.
/// lheap: true if path should be added to heap.
/// lkeep: true if path should be kept for output (only meaningful if lheap).
/// rpath: total path length (computed and returned).
/// xcalcx: running max criterion (passed in/out).
/// iclus[0..natx]: cluster membership (0=inside rfms, 1=outside).
void ccrit(int npat, const int ipat[],
           const float ckspc[], const float fbetac[], const float xlamc[],
           float rmax, float pcrith, float pcritk, int nncrit,
           const int ipot[], float& rpath, bool& lheap, bool& lkeep,
           float& xcalcx, const int iclus[], const AtomData& atoms);

} // namespace feff::path
