// Serial self-consistency loop for muffin-tin potentials.
// Converted from: src/POT/scmt.f
//
// Finds new Fermi level (xmu), electron counts (qnrm), and new valence
// densities (rhoval) through complex-energy contour integration.

#include "scmt.hpp"

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include "../common/logging.hpp"
#include "../common/grid_interpolation.hpp"

#include <cmath>
#include <complex>
#include <cstring>

// Include actual headers for subroutines called by scmt
#include "rholie.hpp"
#include "ff2g.hpp"
#include "coulom.hpp"
#include "../fms/fmsie.hpp"

// Forward declarations for subroutines called by scmt
namespace feff::pot {
    // Defined in grids.cpp
    void grids(double ecv, double xmu, int negx, int& neg,
               std::complex<double>* emg, double* step, int nflrx);
    // Defined in broydn.cpp
    void broydn(int iscmt, double ca1, int nph, const double* xnvmu,
                const int* nr05, const double* xnatph,
                const double* rnrm, double* qnrm,
                double* edenvl, double* rhoval, double* dq);
}

namespace feff::pot {

void scmt(bool verbose, int iscmt, double ecv, int nph, int nat,
          double* vclap, double* edens, double* edenvl,
          double* vtot, double* vvalgs,
          double* rmt, double* rnrm, double* qnrm,
          int ixc, double rhoint, double vint, double& xmu, int jumprm,
          double xnferm, double* xnvmu, double* xnval,
          double x0, double* ri, double dx,
          double* xnatph, double* xion, int iunf, int* iz,
          double* adgc, double* adpc, double* dgc, double* dpc,
          int ihole,
          double* rat, int* iatph, int* iphat,
          int* lmaxsc, double* rhoval, double* xnmues,
          bool& ok,
          double rgrd, int nohole, int nscmt_val, int icoul,
          double ca1, float rfms1, int lfms1)
{
    using Complex = std::complex<double>;
    using CFloat  = std::complex<float>;
    auto& log = feff::common::logger();

    // Local arrays — use static for large arrays to avoid stack overflow.
    // scmt + rholie together need ~2 MB of stack which exceeds the 1 MB default.
    static double dmagx[nrptx], dmag0[251];
    // Fortran: vclap(251,0:nphx), vtot(251,0:nphx), etc.
    static double vtotph[nrptx], vvalph[nrptx], dum[nrptx];
    // dgcn(nrptx,30), dpcn(nrptx,30) work arrays
    static double dgcn[nrptx * 30], dpcn[nrptx * 30];

    // ri05(251), nr05(0:nphx)
    static double ri05[251];
    static int nr05_arr[nphx + 1];
    int* nr05 = nr05_arr;

    // Complex work arrays
    Complex xrhoce[(lx + 1) * (nphx + 1)];
    Complex xrhocp[(lx + 1) * (nphx + 1)];
    Complex xrhole[(lx + 1) * (nphx + 1)];
    Complex yrhoce[251 * (nphx + 1)];
    Complex yrhocp[251 * (nphx + 1)];
    Complex yrhole[251 * (lx + 1) * (nphx + 1)];
    Complex ph[(lx + 1) * (nphx + 1)];

    // Complex energy grid
    constexpr int negx = 80;
    Complex emg[negx];
    constexpr int nflrx = 17;
    double step[nflrx];

    // gtr for FMS
    CFloat gtr[(lx + 1) * (nphx + 1)];

    // dq array for broydn/coulom
    double dq[nphx + 1];

    char slog[512];

    static int ient = 0;
    static feff::fms::FMSData fms_data;  // persistent FMS working data

    bool upok = false;
    int idir = 1;
    ient++;

    // Initialize ri05 grid on first call
    if (ient == 1) {
        for (int i = 0; i < 251; i++) {
            ri05[i] = std::exp(-8.8 + 0.05 * i);
        }
    }

    // nr05[iph] will be overwritten by rholie; this pre-init is just a safety net
    for (int iph = 0; iph <= nph; iph++) {
        nr05[iph] = static_cast<int>((std::log(rnrm[iph]) + 8.8) / 0.05) + 1;
        if (nr05[iph] > 251) nr05[iph] = 251;
        if (nr05[iph] < 1)   nr05[iph] = 1;
    }

    if (verbose) {
        std::snprintf(slog, sizeof(slog),
            "              SCF ITERATION NUMBER%3d  OUT OF%3d", iscmt, nscmt_val);
        log.wlog(slog);
            log.wlog(" Calculating energy and space dependent l-DOS.");
        log.wlog(" It takes time ...");
    }

    // Initialize new valence density
    // rhoval is (251, 0:nphx+1) but we only zero 0:nphx
    for (int iph = 0; iph <= nphx; iph++) {
        for (int ir = 0; ir < 251; ir++) {
            // rhoval[ir + iph * 251]  (0-based)
            rhoval[ir + iph * 251] = 0.0;
        }
    }

    // Get complex energy grid
    int neg = 0;
    grids(ecv, xmu, negx, neg, emg, step, nflrx);

    // Initialize energy loop variables
    int ie = 0;
    double xndifp = 0.0;
    double xndif  = 0.0;
    Complex ee = emg[0];
    Complex ep = static_cast<double>(ee.real());

    // Zero xrhoce and xnmues
    for (int iph = 0; iph <= nphx; iph++) {
        for (int il = 0; il <= lx; il++) {
            // xrhoce[(lx+1)*iph + il]
            xrhoce[il + (lx + 1) * iph] = 0.0;
            // xnmues is (0:lx, 0:nphx), row-major: xnmues[il + (lx+1)*iph]
            xnmues[il + (lx + 1) * iph] = 0.0;
        }
    }
    for (int iph = 0; iph <= nphx; iph++) {
        for (int ir = 0; ir < 251; ir++) {
            yrhoce[ir + 251 * iph] = 0.0;
        }
    }
    int iflr = nflrx;
    int iflrp = nflrx;

    // --- Start energy loop (ie) ---
    // This is the Fortran "25 continue" loop with goto-based control flow.
    // We convert it to a while(true) loop with explicit continue/break.
    Complex eref_val(0.0, 0.0);  // declared outside loop for use after it

energy_loop:
    ie++;

    // Save previous values
    for (int iph = 0; iph <= nph; iph++) {
        for (int il = 0; il <= lx; il++) {
            xrhocp[il + (lx + 1) * iph] = xrhoce[il + (lx + 1) * iph];
        }
        for (int i = 0; i < 251; i++) {
            yrhocp[i + 251 * iph] = yrhoce[i + 251 * iph];
        }
    }

    if (ie == 1 || ie % 20 == 0) {
        if (verbose) {
            std::snprintf(slog, sizeof(slog),
                "     point # %3d  energy = %7.3f", ie, ee.real() * hart);
            log.wlog(slog);
        }
    }

    // Loop over unique potentials
    for (int iph = 0; iph <= nph; iph++) {
        for (int i = 0; i < 251; i++) {
            dmag0[i] = 0.0;
        }

        // fixvar: interpolate edens/vtot to xsect grid
        double vjump = 0.0;
        feff::common::fixvar(rmt[iph], &edens[251 * iph], &vtot[251 * iph],
                             dmag0, vint, rhoint, dx, rgrd, jumprm,
                             vjump, ri, vtotph, dum, dmagx);

        if (ixc % 10 >= 5) {
            int jumprm_tmp = jumprm;
            if (jumprm_tmp > 0) jumprm_tmp = 2;
            feff::common::fixvar(rmt[iph], &edenvl[251 * iph], &vvalgs[251 * iph],
                                 dmag0, vint, rhoint, dx, rgrd, jumprm_tmp,
                                 vjump, ri, vvalph, dum, dmagx);
            // jumprm_tmp reset not needed since it's local
        }

        // fixdsx: interpolate Dirac spinors to xsect grid
        feff::common::fixdsx(iph, dx, rgrd, dgc, dpc, dgcn, dpcn);

        // Compute jri and subtract eref
        int jri = static_cast<int>((std::log(rmt[iph]) + x0) / rgrd) + 2;
        int jri1 = jri + 1;
        // 0-based: vtotph[jri1-1] in Fortran is vtotph[jri1] but arrays are 1-based in Fortran
        // In C++ vtotph is 0-based, jri and jri1 are 1-based Fortran indices
        // jri1 in Fortran = jri+1 (1-based), in C++ = jri (0-based for the +1 offset)
        eref_val = vtotph[jri1 - 1]; // 0-based access for Fortran index jri1
        for (int i = 0; i < jri1; i++) {
            vtotph[i] -= eref_val.real();
        }
        if (ixc >= 5) {
            for (int i = 0; i < jri1; i++) {
                vvalph[i] -= eref_val.real();
            }
        } else {
            for (int i = 0; i < jri1; i++) {
                vvalph[i] = vtotph[i];
            }
        }

        int itmp = 0;
        if (iph == 0 && nohole < 0) itmp = ihole;

        // rholie: compute l-DOS contributions
        // adgc is (10, 30, 0:nphx+1), stride 10*30 per iph
        // xrhole is (0:lx, 0:nphx), stride (lx+1) per iph
        // yrhole is (251, 0:lx, 0:nphx), stride 251*(lx+1) per iph
        // ph is (lx+1, 0:nphx), stride (lx+1) per iph
        rholie(ri05, nr05[iph], rgrd, x0, ri, ee,
               ixc, rmt[iph], rnrm[iph],
               vtotph, vvalph, &xnval[30 * iph], dgcn, dpcn,
               eref_val,
               &adgc[10 * 30 * iph], &adpc[10 * 30 * iph],
               &xrhole[(lx + 1) * iph],
               &xrhoce[(lx + 1) * iph],
               &yrhole[251 * (lx + 1) * iph],
               &yrhoce[251 * iph],
               &ph[(lx + 1) * iph],
               iz[iph], xion[iph], iunf, itmp, lmaxsc[iph]);

    } // end iph loop (label 100)

    // Transform energy variables for FMS
    Complex em_val = static_cast<double>(ee.real());
    Complex eref_fms = static_cast<double>(eref_val.real()) - coni * ee.imag();

    // Initialize gtr
    for (int iph0 = 0; iph0 <= nph; iph0++) {
        for (int il = 0; il <= lx; il++) {
            gtr[il + (lx + 1) * iph0] = CFloat(0.0f, 0.0f);
        }
    }

    // Call FMS if rfms1 > 0
    if (rfms1 > 0.0f) {
        if (lfms1 != 0) {
            int iph0 = 0;
            feff::fms::fmsie(verbose, iph0, nph, lmaxsc, ie, em_val, eref_fms, ph,
                  rfms1, lfms1, nat, iphat, rat, gtr, fms_data);
        } else {
            for (int iph0 = 0; iph0 <= nph; iph0++) {
                feff::fms::fmsie(verbose, iph0, nph, lmaxsc, ie, em_val, eref_fms, ph,
                      rfms1, lfms1, nat, iphat, rat, gtr, fms_data);
            }
        }
    }

    // Compute density and electron counts
    double xntot = 0.0;
    Complex fl(0.0, 0.0), fr(0.0, 0.0);
    for (int iph = 0; iph <= nph; iph++) {
        ff2g(&gtr[(lx + 1) * iph], iph, ie, nr05[iph], xrhoce,
             &xrhole[(lx + 1) * iph], xrhocp,
             ee, ep,
             &yrhole[251 * (lx + 1) * iph],
             &yrhoce[251 * iph],
             &yrhocp[251 * iph],
             &rhoval[251 * iph],
             &xnmues[(lx + 1) * iph], xnatph[iph], xntot,
             iflr, iflrp, fl, fr, iunf);
    }

    if (ie != 1) xndifp = xndif;
    xndif = xntot - xnferm;

    // --- Complex-plane Fermi level search driver ---
    // Decide on next energy point based on floor/direction logic

    Complex eref_val_local; // placeholder for eref used across goto

    if ((ie < neg && ient > 1) || (ient == 1 && ie < nflrx)) {
        ep = ee;
        ee = emg[ie]; // emg is 0-based, ie is already incremented
        if (ie == neg - 1) {
            iflrp = 2;
            iflr  = 1;
        }
        goto energy_loop;
    } else if (ient == 1 && ie == nflrx) {
        upok = false;
        idir = 1;
        if (xntot > xnferm) idir = -1;
        ep = ee;
        ee = ee + static_cast<double>(idir) * step[iflr - 1]; // step is 0-based
        goto energy_loop;
    } else if (ient > 1 && ie == neg) {
        upok = true;
        iflrp = 1;
        iflr  = 1;
        idir = -1;
        if (xntot < xnferm) idir = 1;
        ep = ee;
        ee = ee + static_cast<double>(idir) * step[iflr - 1];
        goto energy_loop;
    } else {
        // Check if the Fermi level is found
        if (iflrp == 1 && iflr == 1 && xndifp * xndif <= 0.0) {
            // Fermi level found
            double xmunew;
            double a;
            if (xndif == 0.0) {
                xmunew = ee.real();
                a = 0.0;
            } else {
                a = xndif / (xndif - xndifp);
                for (int i = 0; i < 4; i++) {
                    Complex fxa_val = a * fl + (1.0 - a) * fr;
                    double bb = Complex((ep - ee) * (fr + fxa_val) / 2.0 +
                                coni * ee.imag() * (fr - fl)).imag();
                    double xndif1 = xndif + a * bb;
                    a = a - xndif1 / bb;
                }
                xmunew = ((1.0 - a) * ee + a * ep).real();
            }

            // Add end-cap corrections to configuration and density
            // Factor 2 for spin degeneracy
            for (int iph = 0; iph <= nph; iph++) {
                for (int il = 0; il <= lx; il++) {
                    if (il <= 2 || iunf != 0) {
                        Complex fl2 = xrhocp[il + (lx + 1) * iph] * 2.0;
                        Complex fr2 = xrhoce[il + (lx + 1) * iph] * 2.0;
                        Complex fxa2 = a * fl2 + (1.0 - a) * fr2;
                        double bb2 = Complex((ep - ee) * (fr2 + fxa2) / 2.0 +
                                    coni * ee.imag() * (fr2 - fl2)).imag();
                        xnmues[il + (lx + 1) * iph] += a * bb2;
                    }
                }
                for (int ir = 0; ir < nr05[iph]; ir++) {
                    Complex fl2 = yrhocp[ir + 251 * iph] * 2.0;
                    Complex fr2 = yrhoce[ir + 251 * iph] * 2.0;
                    Complex fxa2 = a * fl2 + (1.0 - a) * fr2;
                    double bb2 = Complex((ep - ee) * (fr2 + fxa2) / 2.0 +
                                coni * ee.imag() * (fr2 - fl2)).imag();
                    rhoval[ir + 251 * iph] += a * bb2;
                }
            }
            // End of Fermi level found case -- fall through to reporting

            // Report configuration; repeat iteration if found bad counts
            ok = true;
            if (verbose) {
                log.wlog("  Electronic configuration");
                log.wlog("   iph    il      N_el");
            }
            for (int ip = 0; ip <= nph; ip++) {
                for (int il = 0; il <= lx; il++) {
                    if (verbose) {
                        std::snprintf(slog, sizeof(slog), "%6d%6d%9.3f",
                                      ip, il, xnmues[il + (lx + 1) * ip]);
                        log.wlog(slog);
                    }
                    // Check occupation numbers are consistent
                    double diff = std::abs(xnmues[il + (lx + 1) * ip] -
                                          xnvmu[il + (lx + 1) * ip]);
                    if (diff > 13.1 || (il == 2 && diff > 9.1) ||
                        (il == 1 && diff > 5.1) || (il == 0 && diff > 1.95)) {
                        log.wlog(" Found bad counts.");
                        std::snprintf(slog, sizeof(slog),
                            "  Occupation number in getorb is %9.3f",
                            xnvmu[il + (lx + 1) * ip]);
                        log.wlog(slog);
                        log.wlog("  Will repeat this iteration ");
                        // if (ient > 1) ok = false;  // commented in Fortran
                    }
                }
            }

            // If ok, update Fermi level and densities
            if (ok) {
                xmu = xmunew;
                // Find rhoval via Broyden algorithm
                broydn(iscmt, ca1, nph, xnvmu,
                       nr05, xnatph, rnrm, qnrm, edenvl, rhoval, dq);

                // Calculate new vclap - overlap Coulomb potential
                coulom(icoul, nph, nr05, rhoval, edenvl, edens,
                       nat, rat, iatph, iphat, rnrm, dq, iz, vclap);

                // Update edens array
                for (int ip = 0; ip <= nph; ip++) {
                    for (int ir = 0; ir < nr05[ip]; ir++) {
                        edens[ir + 251 * ip] = edens[ir + 251 * ip]
                            - edenvl[ir + 251 * ip] + rhoval[ir + 251 * ip];
                    }
                    for (int ir = nr05[ip]; ir < 251; ir++) {
                        edens[ir + 251 * ip] = 0.0;
                        edenvl[ir + 251 * ip] = 0.0;
                    }
                }
            }
            // Return (Fortran RETURN)
            return;

        } else {
            // Continue search -- goto 25 eventually
            if (iflr == iflrp) {
                // Previous step was horizontal
                if (xndifp * xndif <= 0.0) {
                    // Need to step down
                    upok = false;
                    iflrp = iflr;
                    iflr = iflr - 1;
                    ep = ee;
                    ee = static_cast<double>(ee.real()) +
                         coni * 4.0 * step[iflr - 1];
                } else if (std::abs(xndif) > 10.0 * std::abs(xndif - xndifp) &&
                           upok) {
                    // Need to go up one floor
                    iflrp = iflr;
                    if (iflr < nflrx) {
                        iflr = iflr + 1;
                        ep = ee;
                        ee = static_cast<double>(ee.real()) +
                             coni * 4.0 * step[iflr - 1];
                    } else {
                        ep = ee;
                        ee = ee + static_cast<double>(idir) * step[iflr - 1];
                    }
                } else {
                    // Keep same floor and direction
                    ep = ee;
                    ee = ee + static_cast<double>(idir) * step[iflr - 1];
                }
            } else {
                // Previous step was vertical
                idir = -1;
                if (xndif < 0.0) idir = 1;
                iflrp = iflr;
                ep = ee;
                ee = ee + static_cast<double>(idir) * step[iflr - 1];
            }
            goto energy_loop;
        }
    }
}

} // namespace feff::pot
