// Keep path criterion.
// Converted from: src/PATH/mcritk.f

#include "mcritk.hpp"
#include <cmath>
#include <algorithm>

namespace feff::path {

void mcritk(int npat, const int ipat[], const float ri[], const float beta[],
            const int indbet[], const int ipot[], int nncrit,
            const float fbetac[], const float xlamc[], const float ckspc[],
            float& xout, float& xcalcx) {

    // Not defined if last atom is central atom
    if (ipat[npat - 1] == 0) {
        xout = -1.0f;
        return;
    }

    // Compute output importance factor: sum over p of
    // (product of f(beta)/rho for scatterers) * (cos(beta0)/rho(npat+1))
    // with mean free path factor exp(-rtot/xlam)
    float xcalc = 0.0f;
    float rtot = 0.0f;
    for (int i = 0; i <= npat; ++i) {
        rtot += ri[i];
    }

    for (int icrit = 0; icrit < nncrit; ++icrit) {
        float rho = ri[npat] * ckspc[icrit];
        // Fudge cos(beta0)=0.3 minimum to avoid zero
        float x = std::max(std::abs(beta[npat]), 0.3f) / rho;

        for (int iat = 0; iat < npat; ++iat) {
            rho = ri[iat] * ckspc[icrit];
            int ipot0 = ipot[ipat[iat]];
            x *= fbetac[fbetac_idx(indbet[iat], ipot0, icrit)] / rho;
        }

        x *= std::exp(-rtot / xlamc[icrit]);
        xcalc += x;
    }

    if (xcalcx <= 0.0f) xcalcx = xcalc;
    xout = 100.0f * xcalc / xcalcx;
}

} // namespace feff::path
