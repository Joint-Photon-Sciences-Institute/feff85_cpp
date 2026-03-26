// Prepare scattering amplitude arrays for pathfinder criteria.
// Converted from: src/PATH/prcrit.f

#include "prcrit.hpp"
#include "../common/file_io.hpp"
#include "../math/legendre.hpp"
#include "../par/parallel.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <cmath>
#include <complex>
#include <vector>

namespace feff::path {

void prcrit(int& neout, int& nncrit, int& ik0out,
            float cksp[], float fbeta[],
            float ckspc[], float fbetac[],
            std::string potlb0[],
            float xlam[], float xlamc[]) {

    constexpr double eps = 1.0e-16;

    // Read phase.pad
    auto pd = feff::common::read_xsph("phase.pad");

    int ne = pd.ne1;
    int ik0 = pd.ik0;
    int npot = pd.nph;

    // Extract eref for first spin channel
    std::vector<std::complex<double>> eref(pd.ne);
    for (int ie = 0; ie < pd.ne; ++ie) {
        // eref2(ie, nsp=1) -> eref(ie, 0)
        eref[ie] = pd.eref[ie];  // first spin component
    }

    // Extract phase shifts: ph(ie, il+1, iph) from ph4(ie, -il, 1, iph)
    // pd.ph layout: ph4(ne, 2*ltot+1, nsp, nph+1) flat, Fortran column-major
    // Index: ie + ne * (il+ltot + (2*ltot+1) * (isp + nsp * iph))
    auto ph_idx = [&](int ie, int il_neg, int isp, int iph) -> int {
        return ie + pd.ne * ((il_neg + ltot) + (2*ltot+1) * (isp + pd.nsp * iph));
    };

    // Extract lmax: pd.lmax layout is lmax(ne, 0:nph) flat
    auto lmax_at = [&](int ie, int iph) -> int {
        return pd.lmax[ie + pd.ne * iph];
    };

    // Return outputs
    neout = ne;
    ik0out = ik0;
    for (int i = 0; i <= nphx; ++i) {
        potlb0[i] = pd.potlbl[i];
    }

    // |p| at each energy point, and mean free path
    for (int ie = 0; ie < pd.ne; ++ie) {
        auto cktmp = std::sqrt(2.0 * (pd.em[ie] - eref[ie]));
        cksp[ie] = static_cast<float>(cktmp.real() / feff::bohr);
        // xlam
        xlam[ie] = 1.0e10f;
        if (std::abs(cktmp.imag()) > eps) {
            xlam[ie] = static_cast<float>(1.0 / cktmp.imag());
        }
        xlam[ie] *= static_cast<float>(feff::bohr);
    }

    // Make cos(beta) grid: -40 to 40, 81 points from -1 to 1, spaced 0.025
    double dcosb[2 * nbeta + 1];
    for (int ibeta = -nbeta; ibeta <= nbeta; ++ibeta) {
        dcosb[ibeta + nbeta] = 0.025 * ibeta;
    }
    dcosb[0] = -1.0;           // ibeta = -nbeta
    dcosb[2 * nbeta] = 1.0;   // ibeta = +nbeta

    // Make fbeta for all energy points
    int lmaxp1 = pd.lmaxp1;
    std::vector<double> pl(lmaxp1);

    for (int ibeta = -nbeta; ibeta <= nbeta; ++ibeta) {
        feff::math::cpl0(dcosb[ibeta + nbeta], pl.data(), lmaxp1);

        for (int iii = 0; iii <= npot; ++iii) {
            for (int ie = 0; ie < pd.ne; ++ie) {
                std::complex<double> cfbeta(0.0, 0.0);
                int lm = lmax_at(ie, iii);
                for (int il = 0; il <= lm; ++il) {
                    // ph4(ie, -il, spin=0, iph=iii)
                    auto ph_val = pd.ph[ph_idx(ie, -il, 0, iii)];
                    auto tl = (std::exp(2.0 * feff::coni * ph_val) - 1.0)
                              / (2.0 * feff::coni);
                    cfbeta += tl * pl[il] * static_cast<double>(2 * il + 1);
                }
                fbeta[fbeta_idx(ibeta, iii, ie)] = static_cast<float>(std::abs(cfbeta));
            }
        }
    }

    // Make criterion arrays at selected k-points
    // k = 0, 1, 2, 3, 4, 6, 8, 10, 12 inv-Angstrom
    // ik0 is 0-based (from C++ phmesh/wrxsph)
    int iecrit[necrit];
    iecrit[0] = ik0;
    iecrit[1] = ik0 + 5;
    iecrit[2] = ik0 + 10;
    iecrit[3] = ik0 + 15;
    iecrit[4] = ik0 + 20;
    iecrit[5] = ik0 + 30;
    iecrit[6] = ik0 + 34;
    iecrit[7] = ik0 + 38;
    iecrit[8] = ik0 + 40;

    nncrit = 0;
    for (int ie = 0; ie < necrit; ++ie) {
        if (iecrit[ie] >= pd.ne) break;  // 0-based: valid indices are 0..pd.ne-1
        nncrit = ie + 1;
    }
    if (nncrit == 0) feff::par::par_stop("bad nncrit in prcrit");

    for (int icrit = 0; icrit < nncrit; ++icrit) {
        int ie = iecrit[icrit];  // already 0-based
        ckspc[icrit] = cksp[ie];
        xlamc[icrit] = xlam[ie];
        for (int ibeta = -nbeta; ibeta <= nbeta; ++ibeta) {
            for (int iii = 0; iii <= npot; ++iii) {
                fbetac[fbetac_idx(ibeta, iii, icrit)] = fbeta[fbeta_idx(ibeta, iii, ie)];
            }
        }
    }
}

} // namespace feff::path
