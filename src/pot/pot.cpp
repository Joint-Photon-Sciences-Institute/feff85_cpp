// Main POT computational kernel.
// Converted from: src/POT/pot.f (~680 lines)
//
// Computes potentials for an input atomic cluster, returning data needed
// to compute phase shifts.

#include "pot.hpp"
#include "scmt.hpp"

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include "../common/logging.hpp"
#include "../par/parallel.hpp"
#include "../atom/scfdat.hpp"

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <string>

// Include actual headers for subroutines called by pot
#include "istprm.hpp"
#include "ovrlp.hpp"
#include "afolp.hpp"
#include "corval.hpp"
#include "wpot.hpp"
#include "../atom/utility.hpp"  // for feff::atom::potslw

// Forward declarations for subroutines not yet having separate headers
namespace feff::pot {
    // inipot.cpp
    void inipot(double* dgc, double* dpc, double* edenvl,
                double* vvalgs, double* xnmues);
    // moveh.cpp
    void moveh(int nat, const int* iphat, const int* iz, double* rat);
    // scmtmp.cpp (parallel stub)
    void scmtmp(bool verbose, int npr, int iscmt, double ecv, int nph, int nat,
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
                double ca1, float rfms1, int lfms1,
                std::complex<float>* gtr, std::complex<double>* xrhole,
                std::complex<double>* xrhoce, std::complex<double>* yrhole,
                std::complex<double>* yrhoce);
}

namespace feff::pot {

void pot(bool verbose, double rgrd, int nohole,
         int inters, double totvol, double ecv0, int nscmt, int nmix,
         int ntitle, char title[][80],
         int nat, int nph, int ihole, int iafolp, int ixc,
         int* iphat, double* rat, int* iatph, double* xnatph,
         int* novr, int* iphovr, int* nnovr, double* rovr,
         double* folp0, double* xion, int iunf, int* iz, int ipr1,
         int ispec, int jumprm, int* lmaxsc, int icoul,
         double ca1, float rfms1, int lfms1,
         double& rnrmav, double& xmu, double& vint, double& rhoint,
         double& emu, double& s02, double& erelax, double& wp,
         double& rs, double& xf, double& qtotel,
         int* imt, double* rmt, int* inrm, double* rnrm, double* folpx,
         double* dgc0, double* dpc0,
         double* dgc, double* dpc, double* adgc, double* adpc,
         double* edens, double* vclap, double* vtot,
         double* edenvl, double* vvalgs, double* dmag, double* xnval,
         double* eorb, int* kappa, int* iorb,
         double* qnrm, double* xnmues, int& nhtmp)
{
    using Complex = std::complex<double>;
    using CFloat  = std::complex<float>;
    auto& log = feff::common::logger();
    auto& par = feff::par::state();

    constexpr int Maxprocs = 1;
    constexpr double tolq = 1.0e-3;
    constexpr double tolmu = 3.0e-3;

    char slog[512];

    // Local arrays allocated on the stack
    // rho(251, 0:nphx+1), vcoul(251, 0:nphx+1) - free atom data
    double rho[251 * (nphx + 2)];
    double vcoul[251 * (nphx + 2)];
    double dr[251], drho[251], dvcoul_arr[251];

    // vclapp(251, 0:nphx) - saved Coulomb potential
    double vclapp[251 * (nphx + 1)];

    // rhoval(251, 0:nphx+1)
    double rhoval[251 * (nphx + 2)];

    // folp(0:nphx) - local copy
    double folp[nphx + 1];

    // ri(nrptx) - radial grid
    double ri[nrptx];

    // xnvmu(0:lx, 0:nphx+1) - occupation numbers
    double xnvmu[(lx + 1) * (nphx + 2)];

    // norb(0:nphx+1)
    int norb[nphx + 2];

    // qold(0:nphx)
    double qold[nphx + 1];

    // Work arrays for parallel version (allocated but only used if npr > 1)
    CFloat gtr[(lx + 1) * (nphx + 1) * Maxprocs];
    Complex xrhoce_w[(lx + 1) * (nphx + 1) * Maxprocs];
    Complex xrhole_w[(lx + 1) * (nphx + 1) * Maxprocs];
    Complex yrhoce_w[251 * (nphx + 1) * Maxprocs];
    Complex yrhole_w[251 * (lx + 1) * (nphx + 1) * Maxprocs];

    bool ok = false;
    bool lpass = false;

    // Initialize kappa and eorb arrays
    // Fortran: kappa(30, 0:nphx+1), eorb(30, 0:nphx+1)
    for (int i = 0; i < 30; i++) {
        for (int j = 0; j < nphx + 2; j++) {
            kappa[i + 30 * j] = 0;
            eorb[i + 30 * j] = 0.0;
        }
    }

    // Josh - save nohole value, reset nohole=2 to 0
    nhtmp = nohole;
    if (nohole == 2) nohole = 0;

    // Copy ecv0 and folp0 to local working copies
    double ecv = ecv0;
    for (int i = 0; i <= nph; i++) {
        folp[i] = folp0[i];
    }

    // Initialize arrays
    inipot(dgc, dpc, edenvl, vvalgs, xnmues);

    // Increase length of hydrogen bonds for potential only
    moveh(nat, iphat, iz, rat);

    // Determine if we need two passes (for ION)
    int nfree = 1;
    for (int i = 0; i <= nph; i++) {
        if (std::abs(xion[i]) > 1.0e-3) nfree = 2;
    }

    // === Free atom potentials and densities ===
    // Call twice if any xion != 0 (first time with xion=0 to set rnrm)
    double efrozn = 0.0, etfin = 0.0, etinit = 0.0;

    for (int ifree = 1; ifree <= nfree; ifree++) {
        int ispinr = 0;
        etfin = 0.0;

        for (int iph = 0; iph <= nph; iph++) {
            if (verbose) {
                std::snprintf(slog, sizeof(slog),
                    "    free atom potential and density for atom type%5d", iph);
                log.wlog(slog);
            }

            // Include corehole if absorber (unless nohole)
            int itmp = (iph == 0) ? ihole : 0;

            // Use local variable et to receive total energy from scfdat,
            // matching the Fortran pattern. Only save etfin for the absorber (iph=0).
            double et = 0.0;

            if (nohole >= 0 && iph == 0) {
                // Use nph+1 slot for absorber with core hole
                double xionp = xion[0];
                if (nfree == 2 && ifree == 1) xionp = 0.0;

                int iph_slot = nph + 1;
                feff::atom::scfdat(
                    ipr1, iph_slot, nph, iz[0], itmp, xionp, iunf,
                    &vcoul[251 * iph_slot], &rho[251 * iph_slot],
                    &dmag[251 * iph_slot], &rhoval[251 * iph_slot],
                    ispinr, dgc0, dpc0, dgc, dpc, adgc, adpc,
                    s02, efrozn, et,
                    &xnvmu[(lx + 1) * iph_slot],
                    &xnval[30 * iph_slot],
                    &iorb[8 * iph_slot], norb[iph_slot],
                    &eorb[30 * iph_slot], &kappa[30 * iph_slot]);
            } else {
                double xionp = xion[iph];
                if (nfree == 2 && ifree == 1) xionp = 0.0;

                feff::atom::scfdat(
                    ipr1, iph, nph, iz[iph], itmp, xionp, iunf,
                    &vcoul[251 * iph], &rho[251 * iph],
                    &dmag[251 * iph], &rhoval[251 * iph],
                    ispinr, dgc0, dpc0, dgc, dpc, adgc, adpc,
                    s02, efrozn, et,
                    &xnvmu[(lx + 1) * iph],
                    &xnval[30 * iph],
                    &iorb[8 * iph], norb[iph],
                    &eorb[30 * iph], &kappa[30 * iph]);
            }

            // etfin is absorbing atom final state total energy
            if (iph == 0) etfin = et;
        }

        if (verbose) {
            log.wlog("    initial state energy");
        }

        // Save initial state energy and spinors for core hole orbital
        ispinr = ihole;
        int itmp = 0;

        if (nohole >= 0) {
            int iph = 0;
            double xionp = xion[iph];
            if (nfree == 2 && ifree == 1) xionp = 0.0;

            feff::atom::scfdat(
                ipr1, iph, nph, iz[iph], itmp, xionp, iunf,
                &vcoul[251 * iph], &rho[251 * iph],
                &dmag[251 * iph], &rhoval[251 * iph],
                ispinr, dgc0, dpc0, dgc, dpc, adgc, adpc,
                s02, efrozn, etinit,
                &xnvmu[(lx + 1) * iph],
                &xnval[30 * iph],
                &iorb[8 * iph], norb[iph],
                &eorb[30 * iph], &kappa[30 * iph]);
        } else {
            int iph_slot = nph + 1;
            double xionp = xion[0];
            if (nfree == 2 && ifree == 1) xionp = 0.0;

            feff::atom::scfdat(
                ipr1, iph_slot, nph, iz[0], itmp, xionp, iunf,
                &vcoul[251 * iph_slot], &rho[251 * iph_slot],
                &dmag[251 * iph_slot], &rhoval[251 * iph_slot],
                ispinr, dgc0, dpc0, dgc, dpc, adgc, adpc,
                s02, efrozn, etinit,
                &xnvmu[(lx + 1) * iph_slot],
                &xnval[30 * iph_slot],
                &iorb[8 * iph_slot], norb[iph_slot],
                &eorb[30 * iph_slot], &kappa[30 * iph_slot]);
        }

        // Testing new potential for the final state
        double hx = 0.05;
        double x0_loc = -8.8;
        if (nohole > 0) {
            int idim = 251;
            for (int i = 0; i < idim; i++) {
                dr[i] = std::exp(x0_loc + hx * i);
            }
            if (nohole == 1) {
                for (int i = 0; i < idim; i++) {
                    drho[i] = dgc0[i] * dgc0[i] + dpc0[i] * dpc0[i];
                }
            } else {
                for (int i = 0; i < idim; i++) {
                    // rho(i,0) -> rho[i + 251*0], rhoval(i,0) -> rhoval[i]
                    // rho(i,nph+1) -> rho[i + 251*(nph+1)]
                    drho[i] = dr[i] * dr[i] *
                        (rho[i] - rhoval[i] -
                         rho[i + 251 * (nph + 1)] + rhoval[i + 251 * (nph + 1)]);
                }
            }
            feff::atom::potslw(dvcoul_arr, drho, dr, hx, idim);
            for (int i = 0; i < idim; i++) {
                // Use 1/2 of core-hole as in transition state
                drho[i] = drho[i] / 2.0 / (dr[i] * dr[i]);
            }
        } else {
            for (int i = 0; i < 251; i++) {
                drho[i] = 0.0;
                dvcoul_arr[i] = 0.0;
            }
        }

        // Compute ionization energies
        erelax = -efrozn - (etfin - etinit);
        emu = etfin - etinit;
        fprintf(stderr, "EMU_DEBUG: etfin=%.15e etinit=%.15e emu=%.15e erelax=%.15e efrozn=%.15e s02=%.15e nohole=%d ispinr=%d\n",
                etfin, etinit, emu, erelax, efrozn, s02, nohole, ispinr);

        // Overlap potentials and densities
        for (int iph = 0; iph <= nph; iph++) {
            if (verbose) {
                std::snprintf(slog, sizeof(slog),
                    "    overlapped potential and density for unique potential%5d", iph);
                log.wlog(slog);
            }
            ovrlp(iph, iphat, rat, iatph, novr, iphovr,
                  nnovr, rovr, iz, nat, rho, dmag,
                  rhoval, vcoul, edens, edenvl, vclap, qnrm);
            // Fortran: vclap(1,0) is vclap[0], vcoul(1,0) is vcoul[0]
            if (iph == 0) emu = emu - vclap[0] + vcoul[0];
        }

        if (ifree == 1) {
            // Set the Norman radii
            for (int iph = 0; iph <= nph; iph++) {
                rnrm[iph] = qnrm[iph];
            }
        }

    } // end ifree loop (label 99)

    // Find total charges for istprm
    qtotel = 0.0;
    for (int iph = 0; iph <= nph; iph++) {
        qtotel += (iz[iph] - xion[iph]) * xnatph[iph];
    }

    // Find muffin tin radii, add gsxc to potentials
    if (verbose) {
        log.wlog("    muffin tin radii and interstitial parameters");
    }

    rmt[0] = -1.0;
    xmu = 100.0;
    double xmunew = 0.0;

    if (iafolp >= 0) {
        for (int iph = 0; iph <= nph; iph++) {
            folpx[iph] = folp[iph];
            folp[iph] = 1.0;
        }
    }

    int idmag = 0;
    istprm(nph, nat, iphat, rat, iatph, xnatph,
           novr, iphovr, nnovr, rovr, folp, folpx, iafolp,
           edens, edenvl, idmag,
           dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,
           ixc, rhoint, vint, rs, xf, xmu, xmunew,
           rnrmav, qtotel, inters, totvol);
    xmu = xmunew;

    // Automatic max reasonable overlap
    if (iafolp >= 0) {
        afolp(verbose, nph, nat, iphat, rat, iatph, xnatph,
              novr, iphovr, nnovr, rovr, folp, folpx, iafolp,
              edens, edenvl,
              dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,
              ixc, rhoint, vint, rs, xf, xmu, xmunew,
              rnrmav, qtotel, inters, totvol);
        xmu = xmunew;
    }

    // wp is plasmon frequency in hartrees
    wp = std::sqrt(12.0 * rs / (fa * fa * fa * fa)) * xf * xf / 2.0;

    // Phase shift calculation - Atom r grid
    double dx = 0.05;
    double x0 = 8.8;

    // Initialize charges for SCF
    for (int iph = 0; iph <= nph; iph++) {
        qnrm[iph] = 0.0;
        qold[iph] = 0.0;
    }

    // --- Label 100: corval entry point (may be re-entered if ok==false) ---
corval_entry:
    if (nscmt > 0 || (ispec != 0 && ispec < 4)) {
        corval(verbose, ecv, xnvmu, eorb, norb, xnval,
               kappa, rgrd, nohole,
               nph, edens, edenvl, vtot, vvalgs,
               rmt, rnrm, ixc, rhoint, vint, jumprm,
               x0, ri, dx, xion, iunf, iz,
               adgc, adpc, dgc, dpc, ihole, lmaxsc);
    }

    // Find total number of valence electrons
    double xntot = 0.0;
    for (int iph = 0; iph <= nph; iph++) {
        double xnvmup = 0.0;
        for (int i = 0; i <= lx; i++) {
            // xnvmu(i, iph) -> xnvmu[i + (lx+1)*iph]
            xnvmup += xnvmu[i + (lx + 1) * iph];
        }
        xntot += xnatph[iph] * xnvmup;
    }

    // Update vxcval if nonlocal exchange
    if (ixc % 10 >= 5) {
        istprm(nph, nat, iphat, rat, iatph, xnatph,
               novr, iphovr, nnovr, rovr, folp, folpx, iafolp,
               edens, edenvl, idmag,
               dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,
               ixc, rhoint, vint, rs, xf, xmu, xmunew,
               rnrmav, qtotel, inters, totvol);
        xmunew = xmu;
    }

    if (verbose) {
        std::snprintf(slog, sizeof(slog), "    : mu_old= %9.3f", xmu * hart);
        log.wlog(slog);
    }

    // --- Label 140: nmix decrement point ---
nmix_entry:
    nmix--;

    // Number of processors for parallel execution
    int npr = par.numprocs;

    // === SCF iteration loop (label 200) ===
    for (int iscmt = 1; iscmt <= nscmt; iscmt++) {
        // Save Coulomb potential
        for (int ip = 0; ip <= nph; ip++) {
            for (int ir = 0; ir < 251; ir++) {
                vclapp[ir + 251 * ip] = vclap[ir + 251 * ip];
            }
        }

        if (npr <= 1) {
            scmt(verbose, iscmt, ecv, nph, nat, vclap, edens,
                 edenvl, vtot, vvalgs, rmt, rnrm, qnrm,
                 ixc, rhoint, vint, xmunew, jumprm,
                 xntot, xnvmu, xnval,
                 x0, ri, dx, xnatph, xion, iunf, iz,
                 adgc, adpc, dgc, dpc, ihole,
                 rat, iatph, iphat, lmaxsc, rhoval, xnmues, ok,
                 rgrd, nohole, nscmt, icoul, ca1, rfms1, lfms1);
        } else {
            scmtmp(verbose, npr, iscmt, ecv, nph, nat, vclap, edens,
                   edenvl, vtot, vvalgs, rmt, rnrm, qnrm,
                   ixc, rhoint, vint, xmunew, jumprm,
                   xntot, xnvmu, xnval,
                   x0, ri, dx, xnatph, xion, iunf, iz,
                   adgc, adpc, dgc, dpc, ihole,
                   rat, iatph, iphat, lmaxsc, rhoval, xnmues, ok,
                   rgrd, nohole, nscmt, icoul, ca1, rfms1, lfms1,
                   gtr, xrhole_w, xrhoce_w, yrhole_w, yrhoce_w);
        }

        // If not ok, restart from corval (goto 100)
        if (!ok) goto corval_entry;

        // Test self-consistency
        lpass = true;
        if (iscmt < nscmt && iscmt <= 3) lpass = false;

        if (verbose) {
            std::snprintf(slog, sizeof(slog), " mu_new= %9.3f", xmunew * hart);
            log.wlog(slog);
        }
        if (std::abs(xmunew - xmu) > tolmu) lpass = false;
        xmu = xmunew;

        // Print charge transfer
        if (verbose) log.wlog(" Charge transfer:  iph  charge(iph) ");
        for (int iph = 0; iph <= nph; iph++) {
            if (verbose) {
                std::snprintf(slog, sizeof(slog), "     %3d%9.3f%9.3f",
                              iph, -qnrm[iph] + xion[iph], 0.0);
                log.wlog(slog);
            }
            if (std::abs(qnrm[iph] - qold[iph]) > tolq) lpass = false;
            qold[iph] = qnrm[iph];

            // Check self-consistency of charges
            double sum = -qnrm[iph];
            for (int il = 0; il <= lx; il++) {
                sum += xnmues[il + (lx + 1) * iph] - xnvmu[il + (lx + 1) * iph];
            }
            if (std::abs(sum) > 0.05) lpass = false;
        }

        if (iscmt == nscmt || lpass) {
            // Restore total density from previous iteration
            for (int ip = 0; ip <= nph; ip++) {
                for (int ir = 0; ir < 251; ir++) {
                    edens[ir + 251 * ip] = edens[ir + 251 * ip]
                        - rhoval[ir + 251 * ip] + edenvl[ir + 251 * ip];
                    vclap[ir + 251 * ip] = vclapp[ir + 251 * ip];
                }
                // Remember the reported charge transfer
                qnrm[ip] = -qnrm[ip] + xion[ip];
            }
            // Exit SCF loop (goto 210)
            goto scf_exit;
        } else {
            // Update valence density
            for (int ip = 0; ip <= nph; ip++) {
                for (int ir = 0; ir < 251; ir++) {
                    edenvl[ir + 251 * ip] = rhoval[ir + 251 * ip];
                }
            }
        }

        // Recompute potentials
        istprm(nph, nat, iphat, rat, iatph, xnatph,
               novr, iphovr, nnovr, rovr, folp, folpx, iafolp,
               edens, edenvl, idmag,
               dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,
               ixc, rhoint, vint, rs, xf, xmu, xmunew,
               rnrmav, qtotel, inters, totvol);
        xmunew = xmu;

        if (nmix > 0) goto nmix_entry;
    } // end SCF iteration loop (label 200)

    // Suspicious exit: ran out of iterations

scf_exit: // label 210
    if (par.worker) goto par_barrier_exit; // label 400

    if (nohole > 0) {
        // Testing new final state potential
        for (int j = 0; j < 251; j++) {
            edens[j] -= drho[j]; // edens(j,0)
        }
        for (int j = 0; j < 251; j++) {
            vclap[j] -= dvcoul_arr[j]; // vclap(j,0)
        }

        istprm(nph, nat, iphat, rat, iatph, xnatph,
               novr, iphovr, nnovr, rovr, folp, folpx, iafolp,
               edens, edenvl, idmag,
               dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,
               ixc, rhoint, vint, rs, xf, xmu, xmunew,
               rnrmav, qtotel, inters, totvol);
    }

    // Diagnostic output
    if (ipr1 >= 2) {
        // Convert char[][80] titles to std::string array for wpot
        std::string title_str[10];
        for (int i = 0; i < ntitle && i < 10; i++) {
            title_str[i] = std::string(title[i]);
        }
        wpot(nph, edens, imt, inrm,
             rho, vclap, vcoul, vtot, ntitle, title_str);
    }

par_barrier_exit: // label 400
    // Copy modified folp back to input array so wrpot can access it.
    // In Fortran, wrpot receives the local folp (modified by afolp).
    for (int i = 0; i <= nph; i++) {
        folp0[i] = folp[i];
    }

    feff::par::par_barrier();
}

} // namespace feff::pot
