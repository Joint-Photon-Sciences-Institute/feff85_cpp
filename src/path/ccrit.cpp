// Combined path criterion: heap + keep decisions.
// Converted from: src/PATH/ccrit.f

#include "ccrit.hpp"
#include "mrb.hpp"
#include "mcrith.hpp"
#include "mcritk.hpp"
#include <cmath>
#include <vector>

namespace feff::path {

void ccrit(int npat, const int ipat[],
           const float ckspc[], const float fbetac[], const float xlamc[],
           float rmax, float pcrith, float pcritk, int nncrit,
           const int ipot[], float& rpath, bool& lheap, bool& lkeep,
           float& xcalcx, const int iclus[], const AtomData& atoms) {

    // Compute ri and beta (cos(beta))
    float ri[npatx + 1];
    float beta[npatx + 1];
    mrb(npat, ipat, ri, beta, atoms);

    // Total path length
    rpath = 0.0f;
    for (int i = 0; i <= npat; ++i) {
        rpath += ri[i];
    }

    // Quick reject on path length
    if (rpath > rmax) {
        lheap = false;
        lkeep = false;
        return;
    }

    // If last atom is central atom, put in heap but don't keep
    if (ipat[npat - 1] == 0) {
        lheap = true;
        lkeep = false;
        return;
    }

    // Make index into fbetac array (nearest cos(beta) grid point)
    int indbet[npatx + 1];
    for (int i = 0; i <= npat; ++i) {
        float tmp = std::abs(beta[i]);
        int n = static_cast<int>(tmp / 0.025f);
        float del = tmp - n * 0.025f;
        if (del > 0.0125f) ++n;
        if (beta[i] < 0) n = -n;
        indbet[i] = n;
    }

    // Check heap criterion
    if (pcrith > 0.0f) {
        float xheap;
        mcrith(npat, ipat, ri, indbet, ipot, nncrit, fbetac, ckspc, xheap);

        if (xheap >= 0.0f && xheap < pcrith) {
            lheap = false;
            lkeep = false;
            return;
        }
    }
    lheap = true;

    // Check keep criterion
    if (pcritk <= 0.0f) {
        lkeep = true;
    } else {
        float xout;
        mcritk(npat, ipat, ri, beta, indbet, ipot, nncrit,
               fbetac, xlamc, ckspc, xout, xcalcx);
        lkeep = (xout >= pcritk);
    }

    // If path is entirely inside the FMS cluster, do not keep
    int nclus = 0;
    for (int i = 0; i < npat; ++i) {
        nclus += iclus[ipat[i]];
    }
    if (nclus == 0) lkeep = false;
}

} // namespace feff::path
