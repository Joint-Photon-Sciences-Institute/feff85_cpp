// Map overlapped potential/density to muffin-tin spheres.
// Converted from: src/POT/ovp2mt.f
//
// Uses the LU-decomposed overlap matrix from movrlp to decompose
// overlapped potential at MT boundary into single-site contributions.

#include "ovp2mt.hpp"

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include "../common/logging.hpp"
#include "../par/parallel.hpp"

#include <cmath>
#include <complex>
#include <cstring>

// Forward declarations for LAPACK and utility routines
extern "C" {
    // LAPACK: solve linear system using LU factorization (single-precision complex)
    void cgetrs_(const char* trans, const int* n, const int* nrhs,
                 const std::complex<float>* a, const int* lda,
                 const int* ipiv, std::complex<float>* b, const int* ldb,
                 int* info);
}

// Import math utilities used by ovp2mt
#include "../math/interpolation.hpp"
#include "../math/sommerfeld.hpp"

namespace feff::pot {

// Bring math utilities into pot namespace for use by converted Fortran code
using feff::math::terp;
using feff::math::somm2;

void ovp2mt(int nph, double* vtot, int lrewr, double qtot,
            const double* ri, const double* xnatph, const bool* lnear,
            const int* inrm, const int* imt,
            const double* rnrm, const double* rmt,
            std::complex<float>* cmovp, int* ipiv,
            double& vint, int inters)
{
    using CFloat = std::complex<float>;

    // Extract ipot and irav from inters
    int ipot = inters % 2;
    int irav = (inters - ipot) / 2;

    // Prepare cvovp vector from vtot
    // cvovp has dimension istatx_dim = novp*(nphx+1)+1
    CFloat cvovp[istatx_dim];

    int ncp = 0;
    for (int ip1 = 0; ip1 <= nph; ip1++) {
        for (int i = 1; i <= novp; i++) {
            // Fortran: ix1 = imt(ip1) - novp + i  (1-based grid index)
            int ix1 = imt[ip1] - novp + i;  // 1-based grid index
            int ir_idx = ix1 - 1;            // convert to 0-based for C++ array access
            cvovp[ncp] = CFloat(static_cast<float>(vtot[ir_idx + 251 * ip1]), 0.0f);
            if (lrewr == 2) {
                cvovp[ncp] = CFloat(cvovp[ncp].real() - static_cast<float>(vint), 0.0f);
            }
            ncp++;
        }
    }

    // Compute vtotav for each potential
    double vtotav[nphx + 1];
    for (int ip1 = 0; ip1 <= nph; ip1++) {
        double rav;
        if (irav == 1) {
            rav = (rmt[ip1] + rnrm[ip1]) / 2.0;
        } else if (irav == 0) {
            rav = rnrm[ip1];
        } else {
            // ri[imt[ip1]] is 0-based for Fortran ri(imt(ip1)+1)
            rav = ri[imt[ip1]];  // imt is 1-based, ri is 0-based -> ri[imt-1+1] = ri[imt]
        }
        if (lnear[ip1]) rav = ri[imt[ip1]];

        // terp: interpolate vtot at rav
        // Fortran: terp(ri, vtot(1,ip1), inrm(ip1)+2, 3, rav, vtotav(ip1))
        terp(ri, &vtot[251 * ip1], inrm[ip1] + 2, 3, rav, vtotav[ip1]);
    }

    int istatx = istatx_dim;
    char trans[] = "N";  // NotTransposed
    int nrhs = 1;

    // Find parameters for interstitial potential
    if (lrewr > 0) {
        // Dealing with potentials
        if (lrewr == 1) {
            // Additional equation to find vint
            cvovp[ncp] = CFloat(0.0f, 0.0f);
            double bsum = 0.0;
            // Switch from average equation for vint to local one
            int nphlst = (ipot == 0) ? nph : 0;
            for (int iph = 0; iph <= nphlst; iph++) {
                cvovp[ncp] = CFloat(
                    cvovp[ncp].real() + static_cast<float>(vtotav[iph] * xnatph[iph]),
                    0.0f);
                bsum += xnatph[iph];
            }
            cvovp[ncp] = CFloat(cvovp[ncp].real() / static_cast<float>(bsum), 0.0f);
            ncp++;
        }

        // Solve using LU factorization
        int info = 0;
        cgetrs_(trans, &ncp, &nrhs, cmovp, &istatx, ipiv, cvovp, &istatx, &info);
        if (info < 0) {
            feff::par::par_stop("    *** Error in cgetrf");
        }

        if (lrewr == 1) {
            vint = static_cast<double>(cvovp[ncp - 1].real()) / 100.0;
        }

        // Rewrite vtot
        for (int iph = 0; iph <= nph; iph++) {
            for (int i = 1; i <= novp; i++) {
                // Fortran: index1 = imt(iph)-novp+i (1-based grid), index2 = i+novp*iph (1-based cvovp)
                int index1 = imt[iph] - novp + i;    // 1-based grid index
                int index2 = (i - 1) + novp * iph;   // 0-based into cvovp
                vtot[(index1 - 1) + 251 * iph] =
                    static_cast<double>(cvovp[index2].real()) + vint;
            }
            // Use second order extrapolation
            // Fortran: j = imt(iph)+1
            int j = imt[iph]; // 0-based index for Fortran imt(iph)+1
            double vj;
            terp(ri, &vtot[251 * iph], imt[iph], 2, ri[j], vj);
            vtot[j + 251 * iph] = vj;
            for (int jj = imt[iph] + 1; jj < 251; jj++) {
                vtot[jj + 251 * iph] = vint;
            }
        }
    } else {
        // Dealing with density calculations
        // vint is total charge inside MT spheres

        int info = 0;
        cgetrs_(trans, &ncp, &nrhs, cmovp, &istatx, ipiv, cvovp, &istatx, &info);
        if (info < 0) {
            feff::par::par_stop("    *** Error in cgetrf");
        }

        vint = 0.0;
        double crho[251];
        for (int iph = 0; iph <= nph; iph++) {
            // Fortran loop: i=1, imt(iph)+2 (1-based)
            int np = imt[iph] + 2;  // 1-based count
            for (int i = 1; i <= np; i++) {
                int i0 = i - 1;  // 0-based
                if (i < imt[iph] - novp + 1) {
                    crho[i0] = vtot[i0 + 251 * iph] * ri[i0] * ri[i0];
                } else if (i <= imt[iph]) {
                    // Fortran: ix1 = novp*iph + i - imt(iph) + novp (1-based)
                    int ix1 = novp * iph + (i - 1) - (imt[iph] - 1) + novp - 1;
                    // Simplify: ix1 = novp*iph + i - imt[iph] + novp - 1 (0-based)
                    crho[i0] = static_cast<double>(cvovp[ix1].real()) * ri[i0] * ri[i0];
                } else {
                    // Extrapolate
                    double val;
                    terp(ri, crho, imt[iph], 2, ri[i0], val);
                    crho[i0] = val;
                }
            }

            double cdum = 0.0;
            double dpas = 0.05;
            somm2(ri, crho, dpas, cdum, rmt[iph], 0, np);
            vint += xnatph[iph] * cdum;
        }
        vint = qtot - vint;
    }
}

} // namespace feff::pot
