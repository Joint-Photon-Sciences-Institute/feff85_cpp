#include "legendre.hpp"

namespace feff::math {

void cpl0(double x, double pl0[], int lmaxp1) {
    int lmax = lmaxp1 - 1;

    pl0[0] = 1.0;
    if (lmaxp1 > 1) {
        pl0[1] = x;
    }
    for (int il = 2; il <= lmax; ++il) {
        int l = il - 1;
        pl0[il] = ((2 * l + 1) * x * pl0[il - 1] - l * pl0[il - 2]) / static_cast<double>(il);
    }
}

} // namespace feff::math
