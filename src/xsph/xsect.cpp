// Cross-section calculation (absorption and scattering).
// Converted from src/XSPH/xsect.f (~761 lines)
//
// This is the core cross-section calculation module. It computes:
// - Atomic absorption cross-section xsec(E)
// - Normalization xsnorm(E)
// - Reduced multipole matrix elements rkk(E,k) for EXAFS
//
// For each energy point, it:
// 1. Calculates the self-energy via xcpot
// 2. For each multipole transition (E1, M1, E2):
//    a. Solves the Dirac equation for the final state
//    b. Extracts the phase shift at the muffin-tin boundary
//    c. Normalizes the wavefunction
//    d. Computes radial matrix elements
// 3. Combines with angular coefficients (bmat) to get xsec

#include "xsect.hpp"
#include "radint.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include "../exch/xcpot.hpp"
#include "../fovrg/dfovrg.hpp"
#include "../math/bessel.hpp"
#include "../math/phase_amplitude.hpp"
#include "../math/bcoef.hpp"
#include "../math/sommerfeld.hpp"
#include "../math/interpolation.hpp"
#include "../common/physics_utils.hpp"
#include "../common/logging.hpp"
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>

namespace feff::xsph {

void xsect(int ipr2, double dx, double x0, const double ri[],
           int ne, int ne1, int ik0, const FeffComplex em[], double edge,
           int ihole, double emu, double corr,
           const double dgc0[], const double dpc0[], int jnew,
           int ixc, int lreal, double rmt, double rnrm, double xmu,
           int iPl,
           const double vtot[], const double vvalgs[],
           const double edens[], const double dmag[], const double edenvl[],
           const double dgcn[][30], const double dpcn[][30],
           const double adgc[][30], const double adpc[][30],
           FeffComplex xsec[], double xsnorm[], FeffComplex rkk[][8],
           int iz, double xion, int iunf, const double xnval[],
           int izstd, const int iorb[], int l2lp,
           int ipol, int ispin, int le2, double angks,
           const FeffComplex ptz[3][3]) {

    int kinit, linit;
    feff::common::setkap(ihole, kinit, linit);

    int imt = static_cast<int>((std::log(rmt) + x0) / dx) + 1;
    int jri = imt + 1;
    int jri1 = jri + 1;
    if (jri1 > nrptx) throw std::runtime_error("jri > nrptx in xsect");

    int inrm = static_cast<int>((std::log(rnrm) + x0) / dx) + 1;
    int jnrm = inrm + 1;

    // Normalize initial state: <i|i>
    double xp[nrptx], xq[nrptx];
    for (int i = 0; i < nrptx; i++) {
        xp[i] = dpc0[i] * dpc0[i];
        xq[i] = dgc0[i] * dgc0[i];
    }
    double xinorm = 2.0 * linit + 2.0;
    feff::math::somm(ri, xp, xq, dx, xinorm, 0, jnrm);

    // Set up bcoef
    int kind[8], lind[8];
    constexpr int bmat_total = (2 * lx + 1) * 2 * 8 * (2 * lx + 1) * 2 * 8;
    FeffComplex bmat[bmat_total];
    bool ltrace = true;
    feff::math::bcoef(kinit, ipol, ptz, le2, ltrace, ispin, angks,
                      kind, lind, bmat);
    int isp = 0;
    if (ispin == 1) isp = nspx - 1;

    // Zero rkk
    for (int ie = 0; ie < nex; ie++) {
        for (int k1 = 0; k1 < 8; k1++) rkk[ie][k1] = 0.0;
    }
    FeffComplex phx[8] = {};

    // Workspace
    double vxcrmu[nrptx], vxcimu[nrptx], gsrel[nrptx], vvxcrm[nrptx], vvxcim[nrptx];
    FeffComplex p_arr[nrptx], q_arr[nrptx], pn[nrptx], qn[nrptx];
    FeffComplex pp[nrptx], qp[nrptx], pnp[nrptx], qnp[nrptx];
    FeffComplex xrc[nrptx], xnc[nrptx], xrcold[nrptx], xncold[nrptx];
    FeffComplex v[nrptx], vval[nrptx], fscf[nrptx];
    double bf[3][nrptx];

    // Plasmon pole data
    double WpCorr[MxPole], AmpFac[MxPole];
    for (int i = 0; i < MxPole; i++) WpCorr[i] = -1.0e30;

    if (iPl > 0 && ixc == 0) {
        std::ifstream fin("exc.dat");
        if (fin.is_open()) {
            std::string line;
            int ipole = 0;
            while (ipole < MxPole && std::getline(fin, line)) {
                if (line.empty() || line[0] == '#' || line[0] == '*' ||
                    line[0] == 'c' || line[0] == 'C') continue;
                double wp_val, gam_val, amp_val;
                if (std::sscanf(line.c_str(), "%lf %lf %lf", &wp_val, &gam_val, &amp_val) == 3) {
                    double rs_local = std::pow(3.0 / (4.0 * pi * edens[jri + 1]), third);
                    double wp_calc = std::sqrt(3.0 / (rs_local * rs_local * rs_local));
                    WpCorr[ipole] = (wp_val / hart) / wp_calc;
                    AmpFac[ipole] = amp_val;
                    ipole++;
                }
            }
            if (ipole < MxPole) WpCorr[ipole] = -1.0e30;
        }
    }

    int ifirst = 0;
    int index = ixc;
    feff::atom::FovrgState state;

    // Energy loop
    for (int ie = 0; ie < ne; ie++) {
        FeffComplex eref;
        feff::exch::xcpot(0, ie, index, lreal, ifirst, jri,
                          em[ie], xmu, vtot, vvalgs, edens, dmag, edenvl,
                          eref, v, vval, iPl, WpCorr, AmpFac,
                          vxcrmu, vxcimu, gsrel, vvxcrm, vvxcim, rnrm);

        FeffComplex p2 = em[ie] - eref;
        FeffComplex ck = std::sqrt(2.0 * p2 + (p2 * alphfs) * (p2 * alphfs));
        FeffComplex xkmt = rmt * ck;

        int ncycle = (index % 10 < 5) ? 0 : 3;
        double omega = (std::real(em[ie]) - edge) + emu;
        omega = std::max(omega, 0.001 / hart);

        // Bessel functions for x-ray k-vector
        double xk0 = omega * alphfs;
        int ilast = jnrm + 6;
        if (ilast < jnew) ilast = jnew;
        if (ilast > nrptx) ilast = nrptx;
        for (int i = 0; i < ilast; i++) {
            double temp_val = xk0 * ri[i];
            if (std::abs(temp_val) < 1.0) {
                FeffComplex xirf_temp, dum_temp;
                for (int ll = 0; ll <= 2; ll++) {
                    feff::math::bjnser(FeffComplex(temp_val, 0.0), ll, xirf_temp, dum_temp, 1);
                    bf[ll][i] = std::real(xirf_temp);
                }
            } else {
                double x = temp_val;
                double sinx = std::sin(x), cosx = std::cos(x);
                bf[0][i] = sinx / x;
                bf[1][i] = sinx / (x * x) - cosx / x;
                bf[2][i] = sinx * (3.0 / (x * x * x) - 1.0 / x) - 3.0 * cosx / (x * x);
            }
        }

        xsnorm[ie] = 0.0;
        xsec[ie] = 0.0;
        if (std::real(em[ie]) < -10.0) continue;
        if (std::imag(p2) <= 0.0 && std::real(p2) <= 0.0) continue;

        // For TDLDA: fscf = 1 (no screening)
        for (int i = 0; i < nrptx; i++) fscf[i] = 1.0;

        // Multipole loop: E1 (mult=0), M1 (mult=1), E2 (mult=2)
        for (int mult = 0; mult <= 2; mult++) {
            int kx, ks;
            if (mult == 0) { kx = 1; ks = 2; }
            else { kx = 1; ks = 6; if (mult == 2) kx = 2; }
            if (mult > 0 && mult != le2) continue;

            int ilast2 = jnrm + 6;
            if (ilast2 < jnew) ilast2 = jnew;

            for (int kdif = -kx; kdif <= kx; kdif++) {
                if (omega <= 0.0) continue;
                int ind = kdif + ks - 1; // 0-based index (Fortran: kdif+ks)
                int ikap = kind[ind];
                if (ikap == 0) continue;

                // l2lp selection
                if (l2lp == 1 && ((kinit < 0 && ind >= 2) || (kinit > 0 && ind != 2))) continue;
                if (l2lp == -1 && ((kinit < 0 && ind != 2) || (kinit > 0 && ind >= 2))) continue;

                int iold = 0;
                int ic3 = 0;

                // Solve Dirac equation for final state
                int irr = -1;
                FeffComplex pu, qu, p2_local = p2;
                int jlast = ilast2;
                feff::fovrg::dfovrg(ncycle, ikap, rmt, jlast, jri, p2_local, dx,
                    ri, v, vval, dgcn, dpcn, adgc, adpc,
                    xnval, pu, qu, p_arr, q_arr,
                    iz, ihole, xion, iunf, irr, ic3, state);

                int lfin = lind[ind];
                int ilp = lfin - 1;
                if (ikap < 0) ilp = lfin + 1;
                FeffComplex jl, nl, jlp1, nlp1;
                feff::math::exjlnl(xkmt, lfin, jl, nl);
                feff::math::exjlnl(xkmt, ilp, jlp1, nlp1);
                FeffComplex ph0, temp;
                feff::math::phamp(rmt, pu, qu, ck, jl, nl, jlp1, nlp1, ikap, ph0, temp);

                // Normalize final state
                double sign = (ikap > 0) ? 1.0 : -1.0;
                FeffComplex factor = ck * alphfs;
                factor = sign * factor / (1.0 + std::sqrt(1.0 + factor * factor));
                FeffComplex dum1 = 1.0 / std::sqrt(1.0 + factor * factor);
                FeffComplex xfnorm = 1.0 / temp * dum1;
                for (int i = 0; i < ilast2; i++) {
                    p_arr[i] *= xfnorm;
                    q_arr[i] *= xfnorm;
                }

                // Radial matrix element (E1 includes TDLDA real part)
                for (int i = 0; i < ilast2; i++) {
                    pp[i] = p_arr[i] * std::real(fscf[i]);
                    qp[i] = q_arr[i] * std::real(fscf[i]);
                }
                int ifl = 1;
                if (izstd > 0) ifl = -1;
                FeffComplex xirf = 0.0;
                radint(ifl, mult, bf, kinit, dgc0, dpc0, ikap, pp, qp,
                       pn, qn, ri, dx, ilast2, iold, xrc, xnc, xrcold, xncold, xirf);
                if (ifl < 0) xirf *= xk0;

                // Store matrix element
                if (mult == 0 || le2 == mult) {
                    rkk[ie][ind] = xirf;
                    phx[ind] = ph0;
                }
                xsnorm[ie] += (std::real(xirf) * std::real(xirf) +
                              std::imag(xirf) * std::imag(xirf)) / (2.0 * kx + 1);

                FeffComplex aa = -coni * xirf * xirf;
                int bmat_idx = feff::math::bmat_index(0, isp, ind, 0, isp, ind);
                xsec[ie] -= aa * bmat[bmat_idx];

                // Irregular solution for atomic cross-section
                if (std::imag(em[ie]) > 0.0) {
                    irr = 1;
                    pu = (nl * std::cos(ph0) + jl * std::sin(ph0)) * rmt * dum1;
                    qu = (nlp1 * std::cos(ph0) + jlp1 * std::sin(ph0)) * factor * rmt * dum1;

                    FeffComplex p2_irr = p2;
                    jlast = ilast2;
                    feff::fovrg::dfovrg(ncycle, ikap, rmt, jlast, jri, p2_irr, dx,
                        ri, v, vval, dgcn, dpcn, adgc, adpc,
                        xnval, pu, qu, pn, qn,
                        iz, ihole, xion, iunf, irr, ic3, state);

                    // N- irregular solution
                    temp = std::exp(coni * ph0);
                    for (int i = 0; i < ilast2; i++) {
                        pn[i] = coni * p_arr[i] - temp * pn[i];
                        qn[i] = coni * q_arr[i] - temp * qn[i];
                    }
                } else {
                    for (int i = 0; i < ilast2; i++) { pn[i] = 0.0; qn[i] = 0.0; }
                }

                // Double radial integral for xsec
                for (int i = 0; i < ilast2; i++) {
                    pp[i] = p_arr[i] * std::real(fscf[i]);
                    qp[i] = q_arr[i] * std::real(fscf[i]);
                    pnp[i] = pn[i] * std::real(fscf[i]);
                    qnp[i] = qn[i] * std::real(fscf[i]);
                }
                ifl = 2;
                if (izstd > 0) ifl = -2;
                FeffComplex xirf2 = 0.0;
                radint(ifl, mult, bf, kinit, dgc0, dpc0, ikap, pp, qp,
                       pnp, qnp, ri, dx, ilast2, iold, xrc, xnc, xrcold, xncold, xirf2);
                if (ifl < 0) xirf2 *= xk0 * xk0;

                xsec[ie] -= xirf2 * bmat[bmat_idx];
            }
        }

        // Apply prefactor
        if (omega > 0.0) {
            double prefac = 4.0 * pi * alpinv / omega * bohr * bohr;
            xsnorm[ie] *= prefac * 2.0 * std::abs(ck);
            double xnorm = std::sqrt(xsnorm[ie]);
            xsec[ie] *= prefac * 2.0 * ck;

            // Put sqrt(prefactor) into rkk
            ck = std::sqrt(FeffComplex(prefac, 0.0) * (2.0 * ck));
            if (std::imag(ck) < 0.0) ck = -ck;
            for (int k = 0; k < 8; k++) {
                if (xnorm != 0.0) {
                    rkk[ie][k] *= ck / xnorm * std::exp(coni * phx[k]);
                }
            }
        }
    } // end energy loop
}

} // namespace feff::xsph
