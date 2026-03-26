// Integrate l-resolved DOS from FMS Green's function.
// Converted from src/POT/ff2g.f

#include "ff2g.hpp"
#include <feff/constants.hpp>
#include <cmath>

namespace feff::pot {

void ff2g(const std::complex<float>* gtr, int iph, int ie, int ilast,
          FeffComplex* xrhoce, const FeffComplex* xrhole,
          FeffComplex* xrhocp,
          FeffComplex ee, FeffComplex ep,
          const FeffComplex* yrhole, FeffComplex* yrhoce, FeffComplex* yrhocp,
          double* rhoval, double* xnmues,
          double xnatph, double& xntot, int iflr, int iflrp,
          FeffComplex& fl, FeffComplex& fr, int iunf)
{
    constexpr int s251 = 251;

    // Chi from FMS is contained in gtr (convert single to double precision)
    FeffComplex cchi[lx + 1];
    for (int j = 0; j <= lx; ++j) {
        cchi[j] = FeffComplex(static_cast<double>(gtr[j].real()),
                              static_cast<double>(gtr[j].imag()));
    }

    // xrhoce(il, iph) = xrhoce[iph*(lx+1) + il]
    for (int il = 0; il <= lx; ++il) {
        xrhoce[iph * (lx + 1) + il] += cchi[il] * xrhole[il];
        if (ie == 1) xrhocp[iph * (lx + 1) + il] = xrhoce[iph * (lx + 1) + il];
    }

    FeffComplex del = ee - ep;
    FeffComplex der = del;
    // If iflr=1, add/subtract integral from point to real axis
    // Factor 2 below comes from spin degeneracy
    if (iflr == 1) der = der - coni * 2.0 * std::imag(ee);
    if (iflrp == 1) del = del + coni * 2.0 * std::imag(ep);

    for (int il = 0; il <= lx; ++il) {
        if (il <= 2 || iunf != 0) {
            fl += 2.0 * xrhocp[iph * (lx + 1) + il] * xnatph;
            fr += 2.0 * xrhoce[iph * (lx + 1) + il] * xnatph;
            xnmues[il] += std::imag(xrhoce[iph * (lx + 1) + il] * der +
                                     xrhocp[iph * (lx + 1) + il] * del);
            xntot += xnmues[il] * xnatph;
        }
    }

    // Calculate r-dependent l-DOS for later use
    // yrhole[ir * (lx+1) + il] (251 x (lx+1), column-major)
    for (int il = 0; il <= lx; ++il) {
        if (il <= 2 || iunf != 0) {
            for (int ir = 0; ir < ilast; ++ir) {
                yrhoce[ir] += cchi[il] * yrhole[il * s251 + ir];
                if (ie == 1) yrhocp[ir] = yrhoce[ir];
            }
        }
    }

    for (int ir = 0; ir < ilast; ++ir) {
        rhoval[ir] += std::imag(yrhoce[ir] * der + yrhocp[ir] * del);
    }
}

} // namespace feff::pot
