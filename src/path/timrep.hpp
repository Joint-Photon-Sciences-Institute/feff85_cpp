#pragma once
// Time-reversal and standard ordering of paths.
// Converted from: src/PATH/timrep.f

#include "path_data.hpp"

namespace feff::path {

/// Time-order a path and return it in standard order with hash.
/// npat: number of path atoms.
/// ipat[0..npatx]: atom indices (ipat[npat] set to 0 internally).
/// rx, ry, rz[0..npatx-1]: output path coordinates in standard frame.
/// dhash: output hash key.
/// ipol, ispin: polarization/spin flags.
/// evec[3], xivec[3]: polarization vectors.
/// eels: EELS flag (1 = no time reversal).
void timrep(int npat, int ipat[], float rx[], float ry[], float rz[],
            double& dhash, int ipol, int ispin,
            const double evec[3], const double xivec[3], int eels,
            const AtomData& atoms);

} // namespace feff::path
