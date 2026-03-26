// Spectroscopic density of states (s_z version, no projection).
// Converted from src/XSPH/rholsz.f

#include "rholsz.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include "../fovrg/dfovrg.hpp"
#include "../math/bessel.hpp"
#include "../math/phase_amplitude.hpp"
#include "../math/sommerfeld.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace feff::xsph {

void rholsz(double dx, double x0, const double ri[],
            FeffComplex em, int ixc, double rmt, double rnrm,
            const double vtot[], const double vvalgs[],
            const double xnval[],
            const double dgcn[][30], const double dpcn[][30],
            FeffComplex eref,
            const double adgc[][30], const double adpc[][30],
            FeffComplex* xrhole, FeffComplex* xrhoce,
            FeffComplex ph[], int iz, double xion, int iunf,
            int ihole, int lmaxsc) {

    constexpr int nrx = nrptx;
    int lmax = lmaxsc;
    if (lmax > lx) lmax = lx;
    if (iz <= 4) lmax = 2;
    if (iz <= 2) lmax = 1;

    FeffComplex vtotc[nrptx], vvalc[nrptx];
    for (int i = 0; i < nrptx; i++) {
        vtotc[i] = vtot[i];
        vvalc[i] = vvalgs[i];
    }

    int imt = static_cast<int>((std::log(rmt) + x0) / dx) + 1;
    int jri = imt + 1;
    if (jri > nrptx) throw std::runtime_error("jri > nrptx in rholsz");

    int inrm = static_cast<int>((std::log(rnrm) + x0) / dx) + 1;
    int jnrm = inrm + 1;

    int nr05 = static_cast<int>((std::log(rnrm) + x0) / 0.05) + 5;
    if (nr05 > 251) nr05 = 251;
    int ilast = static_cast<int>(std::round((nr05 - 1) * 0.05 / dx)) + 1;
    if (ilast > nrptx) ilast = nrptx;

    for (int i = 0; i < 64; i++) { xrhole[i] = 0.0; xrhoce[i] = 0.0; }
    for (int i = 0; i <= lx; i++) ph[i] = 0.0;

    FeffComplex p2 = em - eref;
    int ncycle = (ixc % 10 < 5) ? 0 : 3;
    FeffComplex ck = std::sqrt(2.0 * p2 + (p2 * alphfs) * (p2 * alphfs));
    FeffComplex xkmt = rmt * ck;

    FeffComplex pr[nrx][2][2], qr[nrx][2][2], pn[nrx][2][2], qn[nrx][2][2];
    FeffComplex xpc[nrx], phm[2][2];

    feff::atom::FovrgState state;

    for (int lll = 0; lll <= lmax; lll++) {
        for (int jd = 0; jd <= 1; jd++) {
            int ikap = (lll + jd) * (jd == 0 ? 1 : -1);
            if (ikap == 0) continue;
            int ilp = lll + 1;
            if (ikap > 0) ilp = lll - 1;
            int im = jd;

            for (int j = 0; j < 2; j++) {
                int ic3 = j;
                if (lll == 0 && ic3 == 1) continue;

                // Regular solution
                int irr = -1;
                FeffComplex pu, qu, p2_local = p2;
                feff::fovrg::dfovrg(ncycle, ikap, rmt, ilast, jri, p2_local, dx,
                    ri, vtotc, vvalc, dgcn, dpcn, adgc, adpc,
                    xnval, pu, qu, &pn[0][im][j], &qn[0][im][j],
                    iz, ihole, xion, iunf, irr, ic3, state);

                FeffComplex jl, nl, jlp1, nlp1;
                feff::math::exjlnl(xkmt, lll, jl, nl);
                feff::math::exjlnl(xkmt, ilp, jlp1, nlp1);
                FeffComplex phx, temp;
                feff::math::phamp(rmt, pu, qu, ck, jl, nl, jlp1, nlp1, ikap, phx, temp);
                if (lll == 0) ph[0] = phx;
                phm[im][j] = phx;

                FeffComplex xfnorm = 1.0 / temp;
                for (int i = 0; i < ilast; i++) {
                    pr[i][im][j] = pn[i][im][j] * xfnorm;
                    qr[i][im][j] = qn[i][im][j] * xfnorm;
                }

                // Irregular solution
                irr = 1;
                FeffComplex factor_val = ck * alphfs;
                factor_val = factor_val / (1.0 + std::sqrt(1.0 + factor_val * factor_val));
                if (ikap < 0) factor_val = -factor_val;
                qu = (nlp1 * std::cos(phx) + jlp1 * std::sin(phx)) * factor_val * rmt;
                pu = (nl * std::cos(phx) + jl * std::sin(phx)) * rmt;

                feff::fovrg::dfovrg(ncycle, ikap, rmt, ilast, jri, p2_local, dx,
                    ri, vtotc, vvalc, dgcn, dpcn, adgc, adpc,
                    xnval, pu, qu, &pn[0][im][j], &qn[0][im][j],
                    iz, ihole, xion, iunf, irr, ic3, state);

                // N- irregular solution
                temp = std::exp(coni * phx);
                for (int i = 0; i < ilast; i++) {
                    pn[i][im][j] = coni * pr[i][im][j] - temp * pn[i][im][j];
                    qn[i][im][j] = coni * qr[i][im][j] - temp * qn[i][im][j];
                }
            } // j

            FeffComplex factor_val = ck * alphfs;
            factor_val = factor_val / (1.0 + std::sqrt(1.0 + factor_val * factor_val));
            if (ikap < 0) factor_val = -factor_val;
            FeffComplex temp_norm = 2.0 * ck / (1.0 + factor_val * factor_val) / pi;

            int j = 0;
            // Diagonal radial integral
            for (int i = 0; i < ilast; i++) {
                xpc[i] = pr[i][im][j] * pr[i][im][j] + qr[i][im][j] * qr[i][im][j];
            }
            FeffComplex xirf = FeffComplex(lll * 2 + 2, 0.0);
            int i0 = jnrm + 1;
            feff::math::csomm2(ri, xpc, dx, xirf, rnrm, i0);
            int idx = (ikap + 4) + 8 * (ikap + 4);
            xrhole[idx] = xirf * temp_norm * std::exp(coni * (phm[im][j] + phm[im][j]));

            // Central atom contribution
            for (int i = 0; i < ilast; i++) {
                xpc[i] = pn[i][im][j] * pr[i][im][j] + qn[i][im][j] * qr[i][im][j];
                xpc[i] -= coni * (pr[i][im][j] * pr[i][im][j] + qr[i][im][j] * qr[i][im][j]);
            }
            xirf = 1.0;
            feff::math::csomm2(ri, xpc, dx, xirf, rnrm, i0);
            xrhoce[idx] = -xirf * temp_norm;

            // Cross terms
            if (ikap < -1) {
                int k1 = ikap + 2 * lll + 1;
                for (int i = 0; i < ilast; i++) {
                    xpc[i] = pr[i][0][j] * pr[i][1][j] + qr[i][0][j] * qr[i][1][j];
                }
                xirf = FeffComplex(lll * 2 + 2, 0.0);
                i0 = jnrm + 1;
                feff::math::csomm2(ri, xpc, dx, xirf, rnrm, i0);
                int idx1 = (ikap + 4) + 8 * (k1 + 4);
                int idx2 = (k1 + 4) + 8 * (ikap + 4);
                xrhole[idx1] = xirf * temp_norm * std::exp(coni * (phm[0][j] + phm[1][j]));
                xrhole[idx2] = xrhole[idx1];

                // ic3=1 cross terms
                j = 1;
                FeffComplex xpm = std::exp(coni * (phm[0][j] - phm[1][j])) / 2.0;
                FeffComplex xmp = std::exp(coni * (phm[1][j] - phm[0][j])) / 2.0;
                for (int i = 0; i < ilast; i++) {
                    xpc[i] = (pn[i][0][j] * pr[i][1][j] + qn[i][0][j] * qr[i][1][j]) * xmp +
                             (pn[i][1][j] * pr[i][0][j] + qn[i][1][j] * qr[i][0][j]) * xpm;
                    xpc[i] -= coni * (xpm + xmp) *
                              (pr[i][0][j] * pr[i][1][j] + qr[i][0][j] * qr[i][1][j]);
                }
                xirf = 1.0;
                feff::math::csomm2(ri, xpc, dx, xirf, rnrm, i0);
                xrhoce[idx1] = -xirf * temp_norm;
                xrhoce[idx2] = xrhoce[idx1];
            }
        } // jd
    } // lll

    // Phase shifts for l>0
    for (int lll = 1; lll <= lmax; lll++) {
        int im = 0, ikap = -lll - 1, irr = -1, ic3 = 1;
        FeffComplex pu, qu, p2_local = p2;
        feff::fovrg::dfovrg(ncycle, ikap, rmt, ilast, jri, p2_local, dx,
            ri, vtotc, vvalc, dgcn, dpcn, adgc, adpc,
            xnval, pu, qu, &pr[0][im][0], &qr[0][im][0],
            iz, ihole, xion, iunf, irr, ic3, state);

        FeffComplex jl, nl, jlp1, nlp1;
        feff::math::exjlnl(xkmt, lll, jl, nl);
        feff::math::exjlnl(xkmt, lll + 1, jlp1, nlp1);
        FeffComplex phx, temp;
        feff::math::phamp(rmt, pu, qu, ck, jl, nl, jlp1, nlp1, ikap, phx, temp);
        ph[lll] = phx;
    }
}

} // namespace feff::xsph
