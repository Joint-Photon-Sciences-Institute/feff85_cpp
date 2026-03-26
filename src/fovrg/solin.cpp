// Inward solution of the Dirac equation (irregular solution).
// Converted from: src/FOVRG/solin.f
//
// INDEXING: All arrays 0-based. Parameters jri, imax, iwkb are Fortran
// 1-based. Fortran arr(k) = C++ arr[k-1].

#include "solin.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include "../../src/math/bessel.hpp"
#include <cmath>
#include <algorithm>

namespace feff::fovrg {

void flatv(double r1, double r2, FeffComplex p1, FeffComplex q1,
           FeffComplex en, FeffComplex vav, int ikap,
           FeffComplex& p2, FeffComplex& q2);

void solin(FeffComplex en, FeffComplex fl, int kap, int jri, int imax,
           int ic3, const FeffComplex vm[], int iwkb,
           DiracWorkspaceComplex& work, MeshParamsComplex& mesh)
{
    // Exact translation of solin.f
    // npi=6 derivative slots for Milne predictor-corrector
    constexpr int npi = 6;
    constexpr double test_val = 1.0e+5;

    double cl = work.cl;
    double hx = mesh.hx;
    int ndor = mesh.ndor;
    int np = mesh.np;
    double* dr = mesh.dr;  // 0-based, dr[i] = Fortran dr(i+1)

    FeffComplex* gg = work.gg;  // large component
    FeffComplex* ag = work.ag;
    FeffComplex* gp = work.gp;  // small component
    FeffComplex* ap = work.ap;
    FeffComplex* dv = work.dv;  // potential
    FeffComplex* eg = work.eg;  // exchange (large)
    FeffComplex* ep = work.ep;  // exchange (small)

    double ccl = cl + cl;
    int ihard = 0;
    FeffComplex ec = en / cl;

    // Fortran: do i=2,ndor; ag(i)=0; ap(i)=0; enddo
    for (int i = 1; i < ndor; i++) {
        ag[i] = 0.0;
        ap[i] = 0.0;
    }

    // Fortran: vmh = cl * dv(jri+1)
    // dv(jri+1) is 1-based index jri+1 => C++ dv[jri]
    FeffComplex vmh = cl * dv[jri];
    FeffComplex ck = std::sqrt(2.0 * (en - vmh) + std::pow(alphfs * (en - vmh), 2.0));

    int il = std::abs(kap);
    if (kap < 0) il = il - 1;
    int ilp = il - 1;
    if (kap < 0) ilp = il + 1;
    int ilx = il + 1;
    if (ilp > il) ilx = ilp + 1;

    double xsign = -1.0;
    if (kap > 0) xsign = 1.0;

    FeffComplex factor_v = ck * alphfs;
    factor_v = xsign * factor_v / (1.0 + std::sqrt(1.0 + factor_v * factor_v));
    FeffComplex dum1 = 1.0 / std::sqrt(1.0 + factor_v * factor_v);

    // Fortran: iflat = min(jri, iwkb)  — NO modification
    int iflat = std::min(jri, iwkb);

    // Hankel function arrays
    FeffComplex jl_arr[ltot + 2], hl_arr[ltot + 2];
    // Local derivative arrays for Milne (1-based in Fortran: dg(1)..dg(npi))
    // C++ 0-based: dg_loc[0]..dg_loc[npi-1]
    FeffComplex dg_loc[npi], dp_loc[npi];
    for (int k = 0; k < npi; k++) { dg_loc[k] = 0.0; dp_loc[k] = 0.0; }

    // Hankel init: Fortran do i = jri, imax (1-based)
    // C++ 0-based: i from jri-1 to imax-1
    for (int i = jri - 1; i <= imax - 1; i++) {
        // Fortran j = iflat + npi - i (all 1-based)
        // C++ i corresponds to Fortran i_f = i + 1
        // j_f = iflat + npi - i_f = iflat + npi - (i + 1)
        int j_f = iflat + npi - (i + 1);

        FeffComplex xkmt = ck * dr[i];
        feff::math::besjh(xkmt, ilx, jl_arr, hl_arr);
        gg[i] = hl_arr[il] * dr[i] * dum1;
        gp[i] = hl_arr[ilp] * dr[i] * dum1 * factor_v;

        // Fortran: if (j.gt.0) — j is 1-based dg/dp index
        if (j_f >= 1 && j_f <= npi) {
            FeffComplex f_val = (ec - dv[i]) * dr[i];
            FeffComplex g_val = f_val + ccl * dr[i];
            FeffComplex c3_val = FeffComplex(ic3, 0.0) * vm[i] / (g_val * g_val);
            // Fortran: dg(j) = ... => C++ dg_loc[j_f - 1]
            dg_loc[j_f - 1] = -(g_val * gp[i] - FeffComplex(kap, 0.0) * gg[i]);
            dp_loc[j_f - 1] = -(FeffComplex(kap, 0.0) * gp[i] - (f_val - c3_val) * gg[i]);
        }
    }

    // flatv inward loop: Fortran do i = jri-1, iflat, -1 (1-based)
    // C++ 0-based: i from jri-2 down to iflat-1
    // NOTE: if iflat >= jri, this loop doesn't execute (matching Fortran)
    for (int i = jri - 2; i >= iflat - 1; i--) {
        int i_f = i + 1;  // Fortran 1-based index
        int j_f = iflat + npi - i_f;

        FeffComplex eph;
        // Fortran: if (i.eq.iwkb) — i is 1-based
        if (i_f == iwkb) {
            // eph = cl*(3*dv(iwkb+1) - dv(iwkb+2))/2
            // dv(iwkb+1) => C++ dv[iwkb], dv(iwkb+2) => C++ dv[iwkb+1]
            eph = cl * (3.0 * dv[iwkb] - dv[iwkb + 1]) / 2.0;
            // Fortran: if (iwkb.eq.jri-1) eph = cl*(dv(i)+dv(i+1))/2
            // dv(i) => C++ dv[i_f - 1] = dv[i], dv(i+1) => C++ dv[i_f] = dv[i+1]
            if (iwkb == jri - 1) eph = cl * (dv[i] + dv[i + 1]) / 2.0;
        } else {
            eph = cl * (dv[i] + dv[i + 1]) / 2.0;
        }
        if (ic3 > 0) {
            double rav = (dr[i] + dr[i + 1]) / 2.0;
            FeffComplex ec_loc = std::pow(rav, 3.0) *
                std::pow(FeffComplex(ccl, 0.0) + (en - eph) / cl, 2.0);
            eph = eph + FeffComplex(ic3, 0.0) * cl / ec_loc * (vm[i] + vm[i + 1]) / 2.0;
        }
        // flatv(dr(i+1), dr(i), gg(i+1), gp(i+1), en, eph, kap, gg(i), gp(i))
        // dr(i+1) => C++ dr[i+1], dr(i) => C++ dr[i]
        // gg(i+1) => C++ gg[i+1], gg(i) => C++ gg[i]
        // Wait: Fortran i is 1-based. dr(i+1) = C++ dr[(i_f+1)-1] = dr[i_f] = dr[i+1]? No.
        // Fortran dr(i+1) where i is the 1-based loop variable.
        // C++ equivalent: dr[(i_f + 1) - 1] = dr[i_f] = dr[i + 1]
        // Fortran dr(i) = C++ dr[i_f - 1] = dr[i]
        // Same for gg. OK, these are consistent since our loop var i = i_f - 1.
        FeffComplex p2, q2;
        flatv(dr[i + 1], dr[i], gg[i + 1], gp[i + 1], en, eph, kap, p2, q2);
        gg[i] = p2;
        gp[i] = q2;

        if (j_f >= 1 && j_f <= npi) {
            FeffComplex f_val = (ec - dv[i]) * dr[i];
            FeffComplex g_val = f_val + ccl * dr[i];
            FeffComplex c3_val = FeffComplex(ic3, 0.0) * vm[i] / (g_val * g_val);
            dg_loc[j_f - 1] = -(g_val * gp[i] - FeffComplex(kap, 0.0) * gg[i] + ep[i]);
            dp_loc[j_f - 1] = -(FeffComplex(kap, 0.0) * gp[i] - (f_val - c3_val) * gg[i] - eg[i]);
        }
    }

    // Milne predictor-corrector going inward
    double a1 = hx * 3.3;
    double a2 = -hx * 4.2;
    double a3 = hx * 7.8;
    double a4 = hx * 14.0 / 45.0;
    double a5 = hx * 64.0 / 45.0;
    double a6 = hx * 24.0 / 45.0;

    // Fortran: do i = iflat, 2, -1 (1-based)
    // C++ 0-based: i from iflat-1 down to 1
    // In the body, Fortran accesses gg(i+5), gg(i+3), dv(i-1), dr(i-1), etc.
    // With i_f = i + 1: gg(i_f + 5) = C++ gg[i_f + 5 - 1] = gg[i + 5]
    // dv(i_f - 1) = C++ dv[i_f - 2] = dv[i - 1]
    // gg(i_f - 1) (= result) = C++ gg[i_f - 2] = gg[i - 1]
    for (int i = iflat - 1; i >= 1; i--) {
        int nit = 0;

        // Predictor: gg(i_f+5) = gg[i+5]
        FeffComplex acp = gg[i + 5] +
            a1 * (dg_loc[npi - 1] + dg_loc[npi - 5]) +
            a2 * (dg_loc[npi - 2] + dg_loc[npi - 4]) +
            a3 * dg_loc[npi - 3];
        FeffComplex bcp = gp[i + 5] +
            a1 * (dp_loc[npi - 1] + dp_loc[npi - 5]) +
            a2 * (dp_loc[npi - 2] + dp_loc[npi - 4]) +
            a3 * dp_loc[npi - 3];

        // Corrector base: gg(i_f+3) = gg[i+3]
        FeffComplex ac = gg[i + 3] +
            a4 * dg_loc[npi - 4] +
            a5 * (dg_loc[npi - 1] + dg_loc[npi - 3]) +
            a6 * dg_loc[npi - 2];
        FeffComplex bc = gp[i + 3] +
            a4 * dp_loc[npi - 4] +
            a5 * (dp_loc[npi - 1] + dp_loc[npi - 3]) +
            a6 * dp_loc[npi - 2];

        // Shift derivative window: Fortran do j=1,npi-1; dg(j)=dg(j+1)
        for (int jj = 0; jj < npi - 1; jj++) {
            dg_loc[jj] = dg_loc[jj + 1];
            dp_loc[jj] = dp_loc[jj + 1];
        }

        // dv(i_f - 1) = dv[i - 1], dr(i_f - 1) = dr[i - 1]
        FeffComplex f_val = (ec - dv[i - 1]) * dr[i - 1];
        FeffComplex g_val = f_val + ccl * dr[i - 1];
        FeffComplex c3_val = FeffComplex(ic3, 0.0) * vm[i - 1] / (g_val * g_val);

        // Iterate corrector (label 64 in Fortran)
        while (true) {
            // ep(i_f - 1) = ep[i - 1], eg(i_f - 1) = eg[i - 1]
            dg_loc[npi - 1] = -(g_val * bcp - FeffComplex(kap, 0.0) * acp + ep[i - 1]);
            dp_loc[npi - 1] = -(FeffComplex(kap, 0.0) * bcp - (f_val - c3_val) * acp - eg[i - 1]);

            // gg(i_f - 1) = gg[i - 1]
            gg[i - 1] = ac + a4 * dg_loc[npi - 1];
            gp[i - 1] = bc + a4 * dp_loc[npi - 1];

            if (std::abs(test_val * (gg[i - 1] - acp)) > std::abs(gg[i - 1]) ||
                std::abs(test_val * (gp[i - 1] - bcp)) > std::abs(gp[i - 1])) {
                if (nit < 40) {
                    acp = gg[i - 1];
                    bcp = gp[i - 1];
                    nit++;
                    continue;
                } else {
                    ihard++;
                }
            }
            break;
        }
    }

    // Zero beyond imax: Fortran do i=imax+1,np => C++ i from imax to np-1
    for (int i = imax; i < np; i++) {
        gg[i] = 0.0;
        gp[i] = 0.0;
    }

    // Development coefficients: Fortran ag(1) = gg(1)*dr(1)**(-fl)
    ag[0] = gg[0] * std::pow(dr[0], -fl.real());
    ap[0] = gp[0] * std::pow(dr[0], -fl.real());
}

} // namespace feff::fovrg
