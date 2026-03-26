// Output criterion: plane-wave importance factor and pathfinder criteria.
// Converted from: src/PATH/outcrt.f

#include "outcrt.hpp"
#include "mrb.hpp"
#include "mcrith.hpp"
#include "mcritk.hpp"
#include "../math/integration.hpp"
#include <cmath>
#include <algorithm>
#include <vector>

namespace feff::path {

void outcrt(int npat, const int ipat[],
            const float ckspc[], int nncrit,
            const float fbetac[], const float xlamc[],
            int ne, int ik0, const float cksp[],
            const float fbeta[], const float xlam[],
            const int ipot[],
            float& xport, float& xheap, float& xheapr,
            float& xout, float& xcalcx,
            const AtomData& atoms) {

    // Get ri and beta
    float ri[npatx + 1];
    float beta[npatx + 1];
    mrb(npat, ipat, ri, beta, atoms);

    // Make index into fbeta array
    int indbet[npatx + 1];
    for (int i = 0; i <= npat; ++i) {
        float tmp = std::abs(beta[i]);
        int n = static_cast<int>(tmp / 0.025f);
        float del = tmp - n * 0.025f;
        if (del > 0.0125f) ++n;
        if (beta[i] < 0.0f) n = -n;
        indbet[i] = n;
    }

    // Make pw importance factor by integrating over all points above edge
    float rtot = 0.0f;
    for (int i = 0; i <= npat; ++i) rtot += ri[i];

    // ik0 and ne are 0-based: ik0 is the 0-based index of the k=0 point,
    // ne is the total count of energy points (loop up to ne-1).
    std::vector<float> xporti(ne, 0.0f);
    for (int ie = ik0; ie < ne; ++ie) {
        float rho = ri[npat] * cksp[ie];
        float crit = std::max(std::abs(beta[npat]), 0.3f) / rho;
        for (int iat = 0; iat < npat; ++iat) {
            rho = ri[iat] * cksp[ie];
            int ipot0 = ipot[ipat[iat]];
            crit *= fbeta[fbeta_idx(indbet[iat], ipot0, ie)] / rho;
        }
        crit *= std::exp(-rtot / xlam[ie]);
        xporti[ie] = std::abs(crit);
    }

    // Integrate from ik0 to ne-1 using trapezoidal rule
    int nmax = ne - ik0;  // number of points in the range [ik0, ne-1]
    feff::math::strap(&cksp[ik0], &xporti[ik0], nmax, xport);

    // Heap criterion for this path and time-reversed path
    xheap = -1.0f;
    xheapr = -1.0f;
    mcrith(npat, ipat, ri, indbet, ipot, nncrit, fbetac, ckspc, xheap);

    // Prepare time-reversed path arrays
    int nleg = npat + 1;
    float ri0[npatx + 1];
    int indbe0[npatx + 1];
    int ipat0[npatx];

    // Reverse ri
    for (int i = 0; i < nleg; ++i) {
        ri0[i] = ri[nleg - 1 - i];
    }
    // Reverse indbet and ipat
    indbe0[nleg - 1] = indbet[nleg - 1];
    for (int i = 0; i < nleg - 1; ++i) {
        indbe0[i] = indbet[nleg - 1 - i];
    }
    for (int i = 0; i < npat; ++i) {
        ipat0[i] = ipat[nleg - 2 - i];
    }

    mcrith(npat, ipat0, ri0, indbe0, ipot, nncrit, fbetac, ckspc, xheapr);

    // Keep criterion
    mcritk(npat, ipat, ri, beta, indbet, ipot, nncrit,
           fbetac, xlamc, ckspc, xout, xcalcx);
}

} // namespace feff::path
