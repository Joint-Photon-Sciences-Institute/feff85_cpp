// Legendre normalization factors.
// Converted from GENFMT/snlm.f

#include "snlm.hpp"
#include <cmath>
#include <algorithm>

namespace feff::genfmt {

double factst(double flg[211]) {
    // afac = 1/64 works with double precision
    constexpr double afac = 1.0 / 64.0;

    flg[0] = 1.0;
    flg[1] = afac;

    for (int i = 2; i <= 210; ++i) {
        flg[i] = flg[i - 1] * i * afac;
    }

    return afac;
}

void snlm(int lmaxp1, int mmaxp1, double xnlm[ltot + 1][mtot + 1]) {
    double flg[211];
    double afac = factst(flg);

    // Initialize xnlm
    for (int il = 0; il < ltot + 1; ++il) {
        for (int im = 0; im < mtot + 1; ++im) {
            xnlm[il][im] = 0.0;
        }
    }

    for (int il = 0; il < lmaxp1; ++il) {
        int mmxp1 = std::min(mmaxp1, il + 1);
        for (int im = 0; im < mmxp1; ++im) {
            int l = il;
            int m = im;
            double cnlm = (2.0 * l + 1.0) * flg[l - m] / flg[l + m];
            cnlm = std::sqrt(cnlm) * std::pow(afac, m);
            xnlm[il][im] = cnlm;
        }
    }
}

} // namespace feff::genfmt
