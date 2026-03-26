// Path hashing for degeneracy checking.
// Converted from: src/PATH/phash.f
//
// Hashing scheme: Uses ~15 significant digits of double precision.
// With max 9 legs, factor^9 fits in the mantissa range.
// iscale strips trailing digits to avoid roundoff hash collisions.

#include "phash.hpp"
#include <cmath>

namespace feff::path {

void phash(int npat, const int ipat[], const float rx[], const float ry[],
           const float rz[], double& dhash, const AtomData& atoms) {

    constexpr int iscale = 1000;
    constexpr double factor = 16.12345678;
    constexpr double facto2 = 8.57654321;

    dhash = 0.0;

    // Position-based hash (0-based arrays)
    for (int j = 0; j < npat; ++j) {
        double xx = std::pow(factor, j);
        dhash += xx * (std::lround(rx[j] * iscale)
                     + std::lround(ry[j] * iscale) * 0.894375
                     + std::lround(rz[j] * iscale) * 0.573498);
    }

    // Potential-type hash
    for (int j = 0; j < npat; ++j) {
        double xx = std::pow(facto2, j);
        dhash += xx * iscale * atoms.ipot[ipat[j]];
    }

    dhash += npat * 40000000.0;
}

} // namespace feff::path
