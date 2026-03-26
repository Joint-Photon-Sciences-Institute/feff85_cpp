// L-resolved DOS and phase shifts via Dirac equation.
// Converted from src/POT/rholie.f
//
// This is a critical routine. It calls dfovrg twice per l-channel:
//   1) Regular solution (irr=-1): outward integration, normalized to
//      rmt*(jl*cos(delta) - nl*sin(delta))
//   2) Irregular solution (irr=+1): inward integration, combined with
//      regular solution to form the Green's function.

#include "rholie.hpp"
#include "../fovrg/dfovrg.hpp"
#include "../math/bessel.hpp"
#include "../math/phase_amplitude.hpp"
#include "../math/sommerfeld.hpp"
#include "../math/interpolation.hpp"
#include <feff/constants.hpp>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace feff::pot {

void rholie(double* ri05, int& nr05, double dx, double x0,
            double* ri, FeffComplex em,
            int ixc, double rmt, double rnrm,
            double* vtot, double* vvalgs,
            const double* xnval, double* dgcn, double* dpcn,
            FeffComplex eref,
            const double* adgc, const double* adpc,
            FeffComplex* xrhole, FeffComplex* xrhoce,
            FeffComplex* yrhole, FeffComplex* yrhoce,
            FeffComplex* ph,
            int iz, double xion, int iunf, int ihole, int lmaxsc)
{
    constexpr int nrx = nrptx;
    constexpr int s251 = 251;

    // Convert real potentials to complex for dfovrg
    // Use static to avoid stack overflow (these arrays are large)
    static FeffComplex vtotc[nrptx], vvalc[nrptx];
    for (int i = 0; i < nrptx; ++i) {
        vtotc[i] = vtot[i];
        vvalc[i] = vvalgs[i];
    }

    // Set lmax
    int lmax = lmaxsc;
    if (lmax > lx) lmax = lx;
    if (iz <= 4) lmax = 2;
    if (iz <= 2) lmax = 1;

    // Set imt and jri (use general Loucks grid)
    // rmt is between imt and jri
    int imt = static_cast<int>((std::log(rmt) + x0) / dx) + 1;
    int jri = imt + 1;
    if (jri > nrptx) throw std::runtime_error("jri > nrptx in rholie");
    int inrm = static_cast<int>((std::log(rnrm) + x0) / dx) + 1;
    int jnrm = inrm + 1;

    // Set limits for tabulations
    nr05 = static_cast<int>((std::log(rnrm) + x0) / 0.05) + 5;
    if (nr05 > s251) nr05 = s251;
    // ilast is the last integration point
    int ilast = static_cast<int>(std::round((nr05 - 1) * 0.05 / dx)) + 1;
    if (ilast > nrptx) ilast = nrptx;

    // Initialize yrhole and yrhoce
    for (int lll = 0; lll <= lx; ++lll) {
        for (int j = 0; j < s251; ++j) {
            yrhole[lll * s251 + j] = 0.0;
        }
    }
    for (int j = 0; j < s251; ++j) {
        yrhoce[j] = 0.0;
    }

    // p2 is 0.5*(complex momentum)^2 referenced to energy-dependent xc
    FeffComplex p2 = em - eref;
    int ncycle = ((ixc % 10) < 5) ? 0 : 3;
    FeffComplex ck = std::sqrt(2.0 * p2 + (p2 * alphfs) * (p2 * alphfs));
    FeffComplex xkmt = rmt * ck;

    // Work arrays for regular and irregular solutions
    // Use static to avoid stack overflow (these arrays + FovrgState exceed 1 MB)
    static FeffComplex pr[nrx], qr[nrx], pn[nrx], qn[nrx];
    static FeffComplex xpc[nrx];

    // Need a FovrgState for dfovrg
    static feff::atom::FovrgState fovrg_state;

    // Transpose dgcn/dpcn from column-major (fixdsx layout: dgcn[iorb*nrptx+j])
    // to row-major (dfovrg layout: dgcn_t[j*30+iorb]).
    // fixdsx stores: dgcn(j, iorb) = dgcn[iorb * nrptx + j]   (iorb=0..29, j=0..nrptx-1)
    // dfovrg expects: dgcn_t[j][iorb] = dgcn_t[j * 30 + iorb]  (j=0..nrptx-1, iorb=0..29)
    static double dgcn_t[nrx * 30], dpcn_t[nrx * 30];
    for (int iorb = 0; iorb < 30; ++iorb) {
        for (int j = 0; j < nrx; ++j) {
            dgcn_t[j * 30 + iorb] = dgcn[iorb * nrx + j];
            dpcn_t[j * 30 + iorb] = dpcn[iorb * nrx + j];
        }
    }

    // Transpose adgc/adpc from Fortran column-major (adgc(i,j) = adgc[j*10+i])
    // to row-major (dfovrg layout: adgc_t[i][j] = adgc_t[i*30+j]).
    // Fortran: adgc(10, 30), column-major: adgc(i, j) at offset j*10 + i  (i=0..9, j=0..29)
    // dfovrg: adgc_t[i][j] at offset i*30 + j  (i=0..9, j=0..29)
    double adgc_t[10 * 30], adpc_t[10 * 30];
    for (int j = 0; j < 30; ++j) {
        for (int i = 0; i < 10; ++i) {
            adgc_t[i * 30 + j] = adgc[j * 10 + i];
            adpc_t[i * 30 + j] = adpc[j * 10 + i];
        }
    }

    for (int lll = 0; lll <= lx; ++lll) {
        if (lll > lmax) {
            ph[lll] = 0.0;
            xrhoce[lll] = 0.0;
            xrhole[lll] = 0.0;
            for (int i = 0; i < s251; ++i) {
                yrhole[lll * s251 + i] = 0.0;
            }
            continue;
        }

        // ---- Regular solution (irr = -1) ----
        int ikap = -1 - lll;
        int irr = -1;
        int ic3 = (lll == 0) ? 0 : 1;
        FeffComplex p2_local = p2;
        FeffComplex pu, qu;

        feff::fovrg::dfovrg(ncycle, ikap, rmt, ilast, jri, p2_local, dx,
                            ri, vtotc, vvalc,
                            reinterpret_cast<const double(*)[30]>(dgcn_t),
                            reinterpret_cast<const double(*)[30]>(dpcn_t),
                            reinterpret_cast<const double(*)[30]>(adgc_t),
                            reinterpret_cast<const double(*)[30]>(adpc_t),
                            xnval,
                            pu, qu, pn, qn,
                            iz, ihole, xion, iunf, irr, ic3,
                            fovrg_state);

        // Compute Bessel functions at xkmt
        FeffComplex jl, nl, jlp1, nlp1;
        feff::math::exjlnl(xkmt, lll, jl, nl);
        feff::math::exjlnl(xkmt, lll + 1, jlp1, nlp1);

        // Phase shift and amplitude
        FeffComplex phx, temp;
        feff::math::phamp(rmt, pu, qu, ck, jl, nl, jlp1, nlp1, ikap, phx, temp);
        ph[lll] = phx;

        // Normalize final state at rmt to rmt*(jl*cos(delta) - nl*sin(delta))
        FeffComplex xfnorm = 1.0 / temp;

        // Normalize regular solution
        for (int i = 0; i < ilast; ++i) {
            pr[i] = pn[i] * xfnorm;
            qr[i] = qn[i] * xfnorm;
        }

        // ---- Irregular solution (irr = +1) ----
        irr = 1;
        pu = ck * alphfs;
        pu = -pu / (1.0 + std::sqrt(1.0 + pu * pu));

        // Set initial conditions for irregular solution at ilast
        qu = (nlp1 * std::cos(phx) + jlp1 * std::sin(phx)) * pu * rmt;
        pu = (nl * std::cos(phx) + jl * std::sin(phx)) * rmt;

        p2_local = p2;  // reset p2
        feff::fovrg::dfovrg(ncycle, ikap, rmt, ilast, jri, p2_local, dx,
                            ri, vtotc, vvalc,
                            reinterpret_cast<const double(*)[30]>(dgcn_t),
                            reinterpret_cast<const double(*)[30]>(dpcn_t),
                            reinterpret_cast<const double(*)[30]>(adgc_t),
                            reinterpret_cast<const double(*)[30]>(adpc_t),
                            xnval,
                            pu, qu, pn, qn,
                            iz, ihole, xion, iunf, irr, ic3,
                            fovrg_state);

        // Set N- irregular solution: N = i*R - H*exp(i*ph)
        temp = std::exp(coni * phx);

        // Calculate Wronskian
        qu = 2.0 * alpinv * temp * (pn[jri - 1] * qr[jri - 1] - pr[jri - 1] * qn[jri - 1]);
        qu = 1.0 / qu / ck;
        // qu should be close to 1

        for (int i = 0; i < ilast; ++i) {
            pn[i] = coni * pr[i] - temp * pn[i] * qu;
            qn[i] = coni * qr[i] - temp * qn[i] * qu;
        }

        // Combine all constant factors
        // Relativistic correction to normalization and factor 2*lll+1
        pu = ck * alphfs;
        pu = -pu / (1.0 + std::sqrt(1.0 + pu * pu));
        temp = (2.0 * lll + 1.0) / (1.0 + pu * pu) / pi * ck * 2.0;

        // Build xpc for integration: |R|^2 = pr^2 + qr^2
        for (int i = 0; i < ilast; ++i) {
            xpc[i] = pr[i] * pr[i] + qr[i] * qr[i];
        }

        // Interpolate from Loucks grid to 0.05 grid for yrhole
        for (int ir = 0; ir < nr05; ++ir) {
            FeffComplex tempc;
            feff::math::terpc(ri, xpc, ilast, 3, ri05[ir], tempc);
            tempc = tempc * temp;
            yrhole[lll * s251 + ir] = tempc;
        }

        // Integrate xrhole using extended Simpson (csomm2)
        FeffComplex xirf = FeffComplex(lll * 2 + 2, 0.0);
        int i0 = jnrm + 1;
        feff::math::csomm2(ri, xpc, dx, xirf, rnrm, i0);
        xrhole[lll] = xirf * temp;

        // Only central atom contribution needs irregular solution
        for (int i = 0; i < ilast; ++i) {
            xpc[i] = pn[i] * pr[i] - coni * pr[i] * pr[i]
                    + qn[i] * qr[i] - coni * qr[i] * qr[i];
        }

        for (int ir = 0; ir < nr05; ++ir) {
            FeffComplex tempc;
            feff::math::terpc(ri, xpc, ilast, 3, ri05[ir], tempc);
            yrhoce[ir] -= temp * tempc;
        }

        xirf = FeffComplex(1.0, 0.0);
        feff::math::csomm2(ri, xpc, dx, xirf, rnrm, i0);
        xrhoce[lll] = -xirf * temp;
    }
}

} // namespace feff::pot
