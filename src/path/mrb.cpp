// Make ri, beta path parameters for criterion calculations.
// Converted from: src/PATH/mrb.f

#include "mrb.hpp"
#include "../math/distance.hpp"
#include <cmath>

namespace feff::path {

float dotcos(const float rm1[3], const float r[3], const float rp1[3]) {
    constexpr float eps = 1.0e-8f;

    float cosb = 0.0f;
    for (int i = 0; i < 3; ++i) {
        cosb += (r[i] - rm1[i]) * (rp1[i] - r[i]);
    }

    float denom = feff::math::sdist(r, rm1) * feff::math::sdist(rp1, r);
    if (denom > eps) {
        cosb = cosb / denom;
    } else {
        cosb = 0.0f;
    }
    return cosb;
}

void mrb(int npat, const int ipat[], float ri[], float beta[],
         const AtomData& atoms) {
    int nleg = npat + 1;

    // Build local ipat0 with central atom appended (ipat0[nleg-1] = 0)
    // ipat0 is 0-based: ipat0[0..npat-1] = ipat[0..npat-1], ipat0[npat] = 0
    std::vector<int> ipat0(nleg);
    for (int i = 0; i < npat; ++i) {
        ipat0[i] = ipat[i];
    }
    ipat0[nleg - 1] = 0;

    for (int ileg = 0; ileg < nleg; ++ileg) {
        // Work with atom j
        int j = ileg;
        int jm1 = j - 1;
        int jp1 = j + 1;
        // Wrap around
        if (jm1 < 0)     jm1 = nleg - 1;
        if (jp1 >= nleg)  jp1 = 0;

        int jat   = ipat0[j];
        int jm1at = ipat0[jm1];
        int jp1at = ipat0[jp1];

        ri[ileg] = feff::math::sdist(atoms.rat[jat].data(), atoms.rat[jm1at].data());
        beta[ileg] = dotcos(atoms.rat[jm1at].data(), atoms.rat[jat].data(),
                            atoms.rat[jp1at].data());
    }
}

} // namespace feff::path
