// Heap path criterion.
// Converted from: src/PATH/mcrith.f

#include "mcrith.hpp"
#include <cmath>

namespace feff::path {

void mcrith(int npat, const int ipat[], const float ri[], const int indbet[],
            const int ipot[], int nncrit,
            const float fbetac[], const float ckspc[],
            float& xheap) {

    // Not defined for ss, triangles, or paths ending at central atom
    if (ipat[npat - 1] == 0 || npat <= 2) {
        xheap = -1.0f;
        return;
    }

    // Calculate xheap: sum over nncrit of
    // f(beta1)*f(beta2)*...*f(beta npat-2) / (rho1*rho2*...*rho npat-1)
    // Compare to sum(1/p), multiply by 100 for percent.
    xheap = 0.0f;
    float spinv = 0.0f;

    for (int icrit = 0; icrit < nncrit; ++icrit) {
        // ri[npat-2] * ckspc[icrit]^(-(npat-1))
        float x = std::pow(ckspc[icrit], -(npat - 1)) * ri[npat - 2];

        for (int i = 0; i < npat - 2; ++i) {
            int ipot0 = ipot[ipat[i]];
            x *= fbetac[fbetac_idx(indbet[i], ipot0, icrit)] / ri[i];
        }

        spinv += 1.0f / ckspc[icrit];
        xheap += x;
    }

    xheap = 100.0f * xheap / spinv;
}

} // namespace feff::path
