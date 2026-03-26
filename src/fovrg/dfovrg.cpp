// Main driver for the complex-energy Dirac equation solver (FOVRG).
// Converted from: src/FOVRG/dfovrg.f

#include "dfovrg.hpp"
#include "inmuac.hpp"
#include "muatcc.hpp"
#include "wfirdc.hpp"
#include "potex.hpp"
#include "solout.hpp"
#include "solin.hpp"
#include "radial_integrals_c.hpp"
#include "../../src/math/bessel.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace feff::fovrg {

// =========================================================================
// flatv — exact solution for flat potential using spherical Bessel functions
// Converted from: flatv subroutine nested in dfovrg.f
// =========================================================================
void flatv(double r1, double r2, FeffComplex p1, FeffComplex q1,
           FeffComplex en, FeffComplex vav, int ikap,
           FeffComplex& p2, FeffComplex& q2)
{
    // ck = sqrt(2*(en-vav) + (alphfs*(en-vav))^2)
    FeffComplex ck = std::sqrt(2.0 * (en - vav) + (alphfs * (en - vav)) * (alphfs * (en - vav)));

    int isign, lp, lq;
    if (ikap < 0) {
        isign = -1;
        lp = -ikap - 1;
        lq = lp + 1;
    } else {
        isign = 1;
        lp = ikap;
        lq = lp - 1;
    }

    FeffComplex a_ck = ck * alphfs;
    FeffComplex factor = FeffComplex(isign, 0.0) * a_ck / (1.0 + std::sqrt(1.0 + a_ck * a_ck));

    // Bessel function arrays (size ltot+2)
    FeffComplex jl[ltot + 2], nl[ltot + 2];

    // Find a and b such that p1 = r1*(a*jl + b*nl), q1 = factor*r1*(a*jl' + b*nl')
    FeffComplex xkr = ck * r1;
    feff::math::besjn(xkr, jl, nl);

    // Fortran: a = isign*ck*xkr * (p1*nl(lq+1) - q1*nl(lp+1)/factor)
    // Fortran jl/nl are 1-based: jl(1)=l=0, jl(2)=l=1, so jl(lp+1) accesses l=lp.
    // C++ jl/nl are 0-based: jl[0]=l=0, jl[1]=l=1, so jl[lp] accesses l=lp.
    FeffComplex a_coeff = FeffComplex(isign, 0.0) * ck * xkr *
                          (p1 * nl[lq] - q1 * nl[lp] / factor);
    FeffComplex b_coeff = FeffComplex(isign, 0.0) * ck * xkr *
                          (q1 * jl[lp] / factor - p1 * jl[lq]);

    // Get values at r2
    xkr = ck * r2;
    feff::math::besjn(xkr, jl, nl);
    p2 = r2 * (jl[lp] * a_coeff + nl[lp] * b_coeff);
    q2 = r2 * factor * (jl[lq] * a_coeff + nl[lq] * b_coeff);
}

// =========================================================================
// dfovrg — main driver
// =========================================================================
void dfovrg(int ncycle, int ikap, double rmt, int& jlast, int jri,
            FeffComplex& p2, double dx, const double ri[],
            FeffComplex vxc[], FeffComplex vxcval[],
            const double dgcn[][30], const double dpcn[][30],
            const double adgc[][30], const double adpc[][30],
            const double xnval[30],
            FeffComplex& pu, FeffComplex& qu,
            FeffComplex ps[], FeffComplex qs[],
            int iz, int ihole, double xion, int iunf, int irr, int ic3,
            FovrgState& state)
{
    auto& orb = state.orb;
    auto& config = state.config;
    auto& scf = state.scf;
    auto& work = state.work;
    auto& lagrange = state.lagrange;
    auto& mesh = state.mesh;
    auto& error = state.error;

    FeffComplex aps[10], aqs[10];
    FeffComplex vm[nrptx];

    // Initialize data and test parameters
    mesh.ndor = 3;
    double cl = alpinv;
    work.cl = cl;

    if (irr > 0) {
        // For irregular solution
        mesh.ndor = 2;
        aps[0] = pu;
        aqs[0] = qu;
        for (int i = 0; i < jri; i++) {
            work.gg[i] = ps[i];
            work.gp[i] = qs[i];
        }
    }

    // Extend potential beyond jri
    for (int i = jri; i < nrptx; i++) {
        vxc[i] = vxc[jri];
        vxcval[i] = vxc[jri];
    }

    orb.ibgp = 10;
    error.numerr = 0;
    scf.nz = iz;
    mesh.hx = dx;
    mesh.idim = 1 + static_cast<int>(std::round(250 * 0.05 / dx));
    if (mesh.idim > nrptx) mesh.idim = nrptx;
    if (mesh.idim % 2 == 0) mesh.idim = mesh.idim - 1;
    int idim = mesh.idim;

    // WKB switch point
    double aa = 0.5;
    double rwkb = aa / dx / std::sqrt(std::abs(2.0 * p2.real() + (p2.real() / cl) * (p2.real() / cl) +
                  2.0 * p2.imag() + (p2.imag() / cl) * (p2.imag() / cl)));
    // More precise: rwkb = aa / dx / |sqrt(2*p2 + (p2/cl)^2)|
    FeffComplex ck_est = std::sqrt(2.0 * p2 + (p2 / cl) * (p2 / cl));
    rwkb = aa / dx / std::abs(ck_est);
    double x0 = 8.8;
    int iwkb = static_cast<int>((std::log(rwkb) + x0) / dx) + 2;
    if (iwkb > idim) iwkb = idim;
    if (iwkb < 10) iwkb = 10;

    // Copy orbital data into state
    for (int j = 0; j < 30; j++) {
        for (int i = 0; i < 10; i++) {
            orb.bg[i][j] = adgc[i][j];
            orb.bp[i][j] = adpc[i][j];
        }
    }
    for (int j = 0; j < 30; j++) {
        for (int i = 0; i < idim; i++) {
            orb.cg[i][j] = dgcn[i][j];
            orb.cp[i][j] = dpcn[i][j];
        }
    }

    // Initialize orbital configuration
    inmuac(ihole, xion, iunf, ikap, state);
    int norb = scf.norb;
    config.nmax[norb - 1] = jlast;
    if (iwkb >= jlast - 1) iwkb = idim;

    // Calculate initial photoelectron orbital using LDA
    diff(vxc, ri, ikap, cl, mesh.hx, jri, vm);
    // Fortran: do i = jri, nrptx => vm(i) = 0
    // vm(jri) in 1-based = vm[jri-1] in 0-based
    for (int i = jri - 1; i < nrptx; i++) {
        vm[i] = FeffComplex(0.0, 0.0);
    }

    wfirdc(p2, config.kap, config.nmax, vxc, ps, qs, aps, aqs,
           irr, ic3, vm, jri, iwkb, state);

    if (error.numerr != 0) {
        throw std::runtime_error("error in wfirdc");
    }
    if (ncycle == 0) goto label_999;

    // Use only core electrons for exchange term
    for (int i = 0; i < norb - 1; i++) {
        config.xnel[i] = config.xnel[i] - xnval[i];
    }

    // Take vxcval at the origin
    work.av[1] = work.av[1] + (vxcval[0] - vxc[0]) / cl;
    for (int i = 0; i < iwkb; i++) {
        work.dv[i] = vxcval[i] / cl;
    }
    // Keep dv = vxc/cl above iwkb (already set by wfirdc)

    {
        int nter = 0;

        // Angular coefficients
        muatcc(xnval, state);

        // Iteration over cycles
        do {
            nter++;
            int jriwkb = std::min(jri, iwkb);

            // Calculate exchange potential
            potex(ps, qs, aps, aqs, jriwkb, state);

            // Solve Dirac equation
            if (irr < 0) {
                solout(p2, FeffComplex(orb.fl[norb - 1], 0.0), aps[0], aqs[0],
                       ikap, jri, config.nmax[norb - 1], ic3, vm, iwkb,
                       work, mesh);
            } else {
                solin(p2, FeffComplex(orb.fl[norb - 1], 0.0), ikap,
                      jri, config.nmax[norb - 1], ic3, vm, iwkb,
                      work, mesh);
            }

            // Copy result
            config.scc[norb - 1] = 1.0;
            for (int i = 0; i < idim; i++) {
                ps[i] = work.gg[i];
                qs[i] = work.gp[i];
            }
            for (int i = 0; i < mesh.ndor; i++) {
                aps[i] = work.ag[i];
                aqs[i] = work.ap[i];
            }
        } while (nter <= ncycle);
    }

label_999:
    if (error.numerr == 0) {
        if (irr < 0) {
            // Need pu, qu for regular solution at muffin tin.
            // vu is the interstitial potential at the first flat point (0-based jri
            // corresponds to Fortran 1-based jri+1 = jri1).
            // The wavefunction is evaluated at the last integrated grid point
            // (one before the interstitial), which is jri-1 in 0-based indexing
            // (Fortran 1-based jri).
            FeffComplex vu = vxc[jri];
            flatv(ri[jri - 1], rmt, ps[jri - 1], qs[jri - 1], p2, vu, ikap, pu, qu);
            jlast = config.nmax[norb - 1];
        }
    } else {
        throw std::runtime_error("error in dfovrg");
    }
}

} // namespace feff::fovrg
