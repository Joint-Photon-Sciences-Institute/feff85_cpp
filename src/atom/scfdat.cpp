// Main SCF driver for single-configuration Dirac-Fock atom.
// Converted from src/ATOM/scfdat.f
//
// Ankudinov, Zabinsky, Rehr, Comp.Phys. Comm. 98, p.359 (1996).
// Modified Desclaux multi-configuration code, written by A. Ankudinov 1996.
//
// INDEX CONVENTION: All orbital indices are 0-based in C++.
// Fortran orbital j=1..norb maps to C++ j=0..norb-1.
// Functions ortdat() and lagdat() expect 1-based arguments (Fortran convention).

#include "scfdat.hpp"
#include "inmuat.hpp"
#include "muatco.hpp"
#include "wfirdf.hpp"
#include "potrdf.hpp"
#include "vlda.hpp"
#include "soldir.hpp"
#include "ortdat.hpp"
#include "lagdat.hpp"
#include "tabrat.hpp"
#include "etotal.hpp"
#include "fpf0.hpp"
#include "s02at.hpp"
#include "radial_integrals.hpp"
#include "utility.hpp"
#include "../common/orbital_data.hpp"
#include "../common/logging.hpp"
#include "../math/sommerfeld.hpp"
#include "../par/parallel.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <cmath>
#include <cstring>

namespace feff::atom {

void scfdat(int ipr1, int iph, int nph, int iz, int ihole, double xion,
            int iunf, double vcoul[251], double srho[251], double dmag[251],
            double srhovl[251], int ispinr,
            double dgc0[251], double dpc0[251],
            double* dgc, double* dpc,
            double* adgc, double* adpc,
            double& s02, double& efrozn, double& eatom,
            double xntot[], double xnval[30],
            int indorb[], int& norbp, double eorb[30], int kappa_out[30])
{
    // Flat array index helpers for Fortran-layout 3D arrays (all 0-based):
    //   dgc(251, 30, 0:nphx)  →  dgc[i + 251*(j + 30*iph_val)]
    //   adgc(10, 30, 0:nphx)  →  adgc[i + 10*(j + 30*iph_val)]
    auto dgc_idx = [](int i, int j, int iph_val) -> int {
        return i + 251 * (j + 30 * iph_val);
    };
    auto adgc_idx = [](int i, int j, int iph_val) -> int {
        return i + 10 * (j + 30 * iph_val);
    };

    // Create and initialize AtomState (replaces all COMMON blocks)
    AtomState state;

    auto& orb    = state.orb;
    auto& config = state.config;
    auto& scf    = state.scf;
    auto& work   = state.work;
    auto& lag    = state.lagrange;
    auto& nuc    = state.nuclear;
    auto& mesh   = state.mesh;
    auto& ang    = state.angular;
    auto& solver = state.solver;
    auto& error  = state.error;
    auto& inel   = state.inelastic;

    double xmag[30] = {};
    double xnvalp[30] = {};

    // Prepare output file if needed
    bool open_16 = false;
    if (ipr1 >= 3 && iph <= nph) {
        if (feff::par::state().master) {
            open_16 = true;
            char buf[64];
            std::snprintf(buf, sizeof(buf), " free atom %d", iph);
            feff::common::logger().wlog(buf);
        }
    }
    // Fortran immediately resets open_16 after the conditional
    open_16 = false;

    // Initialize data and test parameters
    int jfail = 0;
    orb.ibgp = 10;
    error.numerr = 0;
    scf.nz = iz;
    int iholep = 0;
    int indorb_val = 0;

    inmuat(ihole, xion, iunf, xnval, iholep, xmag, indorb_val, state);

    // NOTE: iholep is 1-based (from getorb). Keep it 1-based since s02at
    // and other functions expect it. Use iholep-1 for 0-based array access.

    // idfock = 1: pure Dirac-Fock (the only mode used)
    int idfock = 1;
    for (int i = 0; i < 30; ++i)
        xnvalp[i] = 0.0;

    int ilast = 0;

    // Calculate initial orbitals using Thomas-Fermi model (ido=1)
    if (error.numerr == 0) {
        double a = -xion - 1.0;
        wfirdf(config.en, a, config.nq, config.kap, config.nmax, 1, state);
    }

    // SCF iteration setup
    int niter = 30;
    int j = 0;       // current orbital (0-based, Fortran j=1)
    int ind = 1;
    int nter = 0;
    int norb = scf.norb;
    int norbsc = scf.norbsc;

    for (int i = 0; i < norb; ++i)
        config.scw[i] = 0.0;

    mesh.test1 = scf.testy / scf.rap[0];
    mesh.test2 = scf.testy / scf.rap[1];
    int netir = std::abs(niter) * norb;

    // Angular coefficients
    muatco(xnvalp, scf, ang, config);
    if (error.numerr != 0) goto label_711;

    // ===== SCF iteration loop (label 101) =====
label_101:
    {
        int iort = 0;
        nter++;

        if (niter < 0) {
            // Orthogonalization by Schmidt procedure (label 104)
            goto label_104;
        }
        goto label_105;

    label_104:
        ortdat(j + 1, state);   // ortdat expects 1-based index

    label_105:
        mesh.method = 1;
        // Calculate Lagrange parameters
        if (lag.nre[j] > 0 && lag.ipl != 0)
            lagdat(j + 1, 1, state);  // lagdat expects 1-based index

        // Calculate electron potential
        potrdf(j, state);

        // Add potential due to xc with valence electrons
        vlda(xnval, srho, srhovl, dmag, ilast, idfock, state);

        double e = config.en[j];
        mesh.np = mesh.idim;

        // Resolution of the Dirac equation
        int ifail = 0;
        double ainf = orb.cg[config.nmax[j] - 1][j];  // cg(nmax(j),j) -> 0-based

        soldir(config.en[j], orb.fl[j], orb.bg[0][j], orb.bp[0][j], ainf,
               config.nq[j], config.kap[j], config.nmax[j], ifail,
               work, mesh, solver, error);

        if (ifail != 0 && jfail == 0) jfail = j + 1;  // track 1-based
        if (jfail == j + 1 && ifail == 0) jfail = 0;

        if (error.numerr != 0) {
            if (iort != 0 || niter < 0) goto label_711;
            iort = 1;
            goto label_104;
        }

        // label 111: successful Dirac solution
        config.sce[j] = std::abs((e - config.en[j]) / config.en[j]);

        // Find maximum variation of the wave function
        int k = config.nmax[j];
        double pr = 0.0;
        double a_coef = 0.0, b_coef = 0.0;
        for (int i = 0; i < k; ++i) {
            double w = orb.cg[i][j] - work.dg[i];
            if (std::abs(w) > std::abs(pr)) {
                pr = w;
                a_coef = orb.cg[i][j];
                b_coef = work.dg[i];
            }
            w = orb.cp[i][j] - work.dp[i];
            if (std::abs(w) > std::abs(pr)) {
                pr = w;
                a_coef = orb.cp[i][j];
                b_coef = work.dp[i];
            }
        }

        // Acceleration of convergence
        b_coef = config.scc[j];
        cofcon(a_coef, b_coef, pr, config.scw[j]);
        config.scc[j] = b_coef;

        for (int i = 0; i < k; ++i) {
            work.dg[i] = b_coef * work.dg[i] + a_coef * orb.cg[i][j];
            work.dp[i] = b_coef * work.dp[i] + a_coef * orb.cp[i][j];
        }
        for (int i = 0; i < mesh.ndor; ++i) {
            work.ag[i] = b_coef * work.ag[i] + a_coef * orb.bg[i][j];
            work.ap[i] = b_coef * work.ap[i] + a_coef * orb.bp[i][j];
        }

        // Normalization of the wave function
        // Fortran: a = dsordf(j, k, 0, 4, fl(j)) where j is 1-based orbital, k=nmax(j).
        // When jnd=4, dsordf treats its 2nd argument as a POINT COUNT (max0=j_arg),
        // NOT an array index. Since nmax is a 1-based count from Fortran soldir,
        // we pass k directly — the C++ dsordf loop "for l=0; l<max0" covers
        // max0 points (0-based indices 0..max0-1), matching Fortran's "do l=1,max0".
        // BUG FIX: Previously passed k-1, which integrated one point too few,
        // causing incorrect normalization and wrong SCF eigenvalues.
        double a_norm = dsordf(j, k, 0, 4, orb.fl[j], orb, work, config, mesh);
        a_norm = std::sqrt(a_norm);

        for (int i = 0; i < mesh.np; ++i) {
            orb.cg[i][j] = work.dg[i] / a_norm;
            orb.cp[i][j] = work.dp[i] / a_norm;
        }
        for (int i = 0; i < mesh.ndor; ++i) {
            orb.bg[i][j] = work.ag[i] / a_norm;
            orb.bp[i][j] = work.ap[i] / a_norm;
        }

        // Determination of the next orbital to calculate
        // Fortran: if (nter.lt.norbsc .or. (ind.lt.0 .and. j.lt.norbsc))
        // where j is 1-based. In 0-based: j < norbsc-1
        if (nter < norbsc || (ind < 0 && j < norbsc - 1)) {
            j++;
            goto label_451;
        }

        // Fortran: j = j + 1 (advance past current)
        j++;

        // Find orbital with largest wave function error
        pr = 0.0;
        for (int i = 0; i < norbsc; ++i) {
            double w = std::abs(config.scw[i]);
            if (w > pr) {
                pr = w;
                j = i;
            }
        }
        // Fortran: if (j .gt. norbsc) j = 1  → 0-based: if (j >= norbsc) j = 0
        if (j >= norbsc) j = 0;
        if (pr > scf.testy) goto label_421;

        // Find orbital with largest energy error
        pr = 0.0;
        for (int i = 0; i < norbsc; ++i) {
            double w = std::abs(config.sce[i]);
            if (w > pr) {
                pr = w;
                j = i;
            }
        }
        if (pr >= scf.teste) goto label_421;

        // Both converged
        if (ind < 0) goto label_999;
        ind = -1;
        j = 0;
        goto label_451;
    }

label_421:
    ind = 1;

label_451:
    if (nter <= netir) goto label_101;

    // Number of iterations exceeded the limit
    error.numerr = 192011;
    std::strncpy(error.dlabpr, "  scfdat", 8);

label_711:
    messer(error);
    feff::par::par_stop("SCFDAT-1");

label_999:
    norb = scf.norb;   // refresh after potential changes

    if (error.numerr == 0) {
        if (jfail != 0) {
            feff::par::par_stop(
                "  Failed to match lower component, results are meaningless");
        }

        // Tabulation of results
        if (ipr1 >= 5 && iph <= nph) tabrat(state);

        // Total energy
        etotal(16, eatom, state);

        // Prepare for XC energy contribution
        for (int ix = 0; ix < 251; ++ix)
            dmag[ix] = 0.0;
        ilast = 1;

        vlda(xnval, srho, srhovl, dmag, ilast, idfock, state);

        double ecorr = 2.0;
        feff::math::somm(mesh.dr, dmag, dmag, mesh.hx, ecorr, 0, mesh.idim);
        eatom = eatom - ecorr / 4.0;

        // Prepare information for SCMT and core-valence separation
        norbp = norb;
        for (int i = 0; i <= feff::lx; ++i)
            xntot[i] = 0.0;

        for (int jj = 0; jj < norb; ++jj) {
            eorb[jj] = config.en[jj];
            kappa_out[jj] = config.kap[jj];
            int ii = config.kap[jj];
            if (config.kap[jj] < 0) ii = -config.kap[jj] - 1;
            if (ii <= feff::lx) xntot[ii] += xnval[jj];
        }

        // Get difference in spin-up and -down densities per spin
        double spin = 0.0;
        for (int i = 0; i < mesh.idim; ++i)
            dmag[i] = 0.0;

        for (int iorb = 0; iorb < norb; ++iorb) {
            spin += xmag[iorb];
            for (int i = 0; i < mesh.np; ++i) {
                dmag[i] += xmag[iorb] *
                    (orb.cg[i][iorb] * orb.cg[i][iorb] +
                     orb.cp[i][iorb] * orb.cp[i][iorb]);
            }
        }
        if (spin > 0.0) {
            for (int i = 0; i < mesh.np; ++i)
                dmag[i] /= spin;
        }

        // Return Coulomb potential
        potslw(vcoul, srho, mesh.dr, mesh.hx, mesh.idim);
        for (int i = 0; i < 251; ++i)
            vcoul[i] = vcoul[i] - static_cast<double>(scf.nz) / mesh.dr[i];

        // Return srho as 4*pi*density instead of 4*pi*density*r**2
        for (int i = 0; i < 251; ++i) {
            double r2 = mesh.dr[i] * mesh.dr[i];
            srho[i]   /= r2;
            dmag[i]   /= r2;
            srhovl[i] /= r2;
        }

        if (ipr1 >= 3 && iph <= nph) {
            // Close output file (in Fortran this closes unit 16)
        }

        // Save central atom Dirac components for spin-polarized case
        if (ispinr != 0) {
            // Need kap(i) for central atom without core hole; all output of
            // getorb is dummy except iholep and kap(i) (stored in nq(i))
            int dummy_norb, dummy_norbco, dummy_iorb;
            int dummy_nqn[30], dummy_nk[30];
            double dummy_xnel[30], dummy_xnval[30], dummy_xmag[30];

            feff::common::getorb(iz, ispinr, xion, iunf,
                                 dummy_norb, dummy_norbco, dummy_iorb,
                                 iholep, dummy_nqn, dummy_nk,
                                 dummy_xnel, dummy_xnval, dummy_xmag);

            // iholep is 1-based from getorb. Keep as-is for s02at which expects 1-based.
            // Use (iholep-1) for 0-based array access.

            // Store no-hole kappa values in nq for s02 overlap comparison
            // (Fortran puts these into COMMON nq via the getorb call)
            for (int i = 0; i < 30; ++i)
                config.nq[i] = dummy_nk[i];

            // Copy central atom wave functions for core-hole orbital
            // iholep is 1-based; use iholep-1 for 0-based C++ arrays
            int ihp0 = iholep - 1;  // 0-based index
            for (int i = 0; i < config.nmax[ihp0]; ++i) {
                dgc0[i] = orb.cg[i][ihp0];
                dpc0[i] = orb.cp[i][ihp0];
            }
            for (int i = config.nmax[ihp0]; i < 251; ++i) {
                dgc0[i] = 0.0;
                dpc0[i] = 0.0;
            }
        }

        // Copy all orbital wave functions to output arrays
        for (int jj = 0; jj < 30; ++jj) {
            int nm = (jj < norb) ? config.nmax[jj] : 0;
            for (int i = 0; i < nm; ++i) {
                dgc[dgc_idx(i, jj, iph)] = orb.cg[i][jj];
                dpc[dgc_idx(i, jj, iph)] = orb.cp[i][jj];
            }
            for (int i = nm; i < 251; ++i) {
                dgc[dgc_idx(i, jj, iph)] = 0.0;
                dpc[dgc_idx(i, jj, iph)] = 0.0;
            }
            for (int i = 0; i < 10; ++i) {
                adgc[adgc_idx(i, jj, iph)] = orb.bg[i][jj];
                adpc[adgc_idx(i, jj, iph)] = orb.bp[i][jj];
            }
        }
    }

    // Calculate overlap integrals for final and initial state orbitals
    // of the central atom. This runs only in the last call of scfdat
    // from subroutine pot (ihole=0 and iholep=ispinr != 0).
    if (iholep > 0 && iholep < 30 && ihole <= 0) {
        efrozn = config.en[iholep - 1];  // iholep is 1-based; arrays are 0-based
        double ovpint[30][30] = {};

        for (int i = 0; i < norb; ++i) {
            // Handle special case when electron added to new orbital
            // After getorb call above, nq holds kappa for non-hole atom
            int itr;
            if (config.nq[i] == config.kap[i]) {
                itr = 0;
            } else if (i + 1 < 30 && config.nq[i + 1] == config.kap[i]) {
                itr = 1;
            } else {
                feff::common::logger().wlog(
                    "  If it is not la, gd or np, please, give us a call");
                feff::common::logger().wlog("  s02 is overestimated");
                for (int jj = 0; jj < i; ++jj)
                    ovpint[jj][i] = 0.0;
                ovpint[i][i] = 1.0;
                goto label_780;
            }

            {
                int i0 = i + itr;  // 0-based
                int iph1 = 0;
                if (iph == 0) iph1 = nph + 1;

                // Load initial-state orbital into Dirac workspace
                for (int ir = 0; ir < mesh.idim; ++ir) {
                    work.dg[ir] = dgc[dgc_idx(ir, i0, iph1)];
                    work.dp[ir] = dpc[dgc_idx(ir, i0, iph1)];
                }
                for (int ir = 0; ir < mesh.ndor; ++ir) {
                    work.ag[ir] = adgc[adgc_idx(ir, i0, iph1)];
                    work.ap[ir] = adpc[adgc_idx(ir, i0, iph1)];
                }

                // Compute overlap integrals with all same-kappa orbitals
                for (int jj = 0; jj < norb; ++jj) {
                    if (config.kap[i] != config.kap[jj]) continue;
                    // dsordf with jnd=3: overlap of orbital jj with workspace
                    ovpint[i][jj] = dsordf(jj, jj, 0, 3, orb.fl[i],
                                           orb, work, config, mesh);
                }
            }

        label_780:;
        }

        // Subtract valence occupations to get core-only
        for (int jj = 0; jj < norb; ++jj)
            config.xnel[jj] -= xnval[jj];

        // fpf0 call is commented out in the Fortran source

        // Calculate S02
        s02at(iholep, norb, config.kap, config.xnel, ovpint, s02);
    }

}

} // namespace feff::atom
