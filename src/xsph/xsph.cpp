// Master orchestrator for phase shifts and cross-sections.
// Converted from src/XSPH/xsph.f (~570 lines)

#include "xsph.hpp"
#include "getedg.hpp"
#include "phmesh.hpp"
#include "phase.hpp"
#include "xsect.hpp"
#include "wrxsph.hpp"
#include "axafs.hpp"
#include "wphase.hpp"
#include "szlz.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <feff/json_io.hpp>
#include "../common/grid_interpolation.hpp"
#include "../common/physics_utils.hpp"
#include "../common/logging.hpp"
#include "../pot/istprm.hpp"
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <fstream>
#include <vector>

namespace feff::xsph {

void xsph(bool wrxsec, bool verbose, const std::string& phpad,
           int ipr2, int ispec, double vixan, double xkstep, double xkmax,
           double gamach, double rgrd,
           int nph, int lmaxph[], const char potlbl[][7],
           const double spinph[], const int iatph[], int nat,
           const double* rat, const int iphat[],
           int ixc, double vr0, double vi0, int ixc0, int lreal,
           float rfms2, int lfms2, int l2lp,
           int ipol, int ispin, int le2, double angks, const FeffComplex ptz[3][3],
           int iPl, int izstd, int ifxc, int ipmbse, int itdlda, int nonlocal,
           int ntitle, const char title[][80], double rnrmav,
           double xmu, double vint, double rhoint,
           double emu, double s02, double erelax, double wp, double ecv,
           double rs, double xf, double qtotel,
           const int imt[], const double rmt[], const int inrm[],
           const double rnrm[], const double folp[], double folpx[],
           const double xnatph[],
           const double dgc0_in[], const double dpc0_in[],
           const double* dgc, const double* dpc,
           const double* adgc, const double* adpc,
           double* edens, double* vclap, double* vtot,
           double* edenvl, double* vvalgs, double* dmag,
           const double* xnval, const int* iorb,
           int nohole, int ihole,
           int inters, double totvol, int iafolp,
           const double xion[], int iunf, const int iz[], int jumprm) {

    // Phase shift calculation grids
    double dx = 0.05;
    double x0 = 8.8;
    double dxnew = rgrd;

    // Optical constants edge energy lookup (disabled by default)
    bool lopt = false;
    int ik0 = 0;
    if (lopt) {
        getedg(ihole, iz[0], emu);
        ik0 = 0;
        if (verbose) {
            char slog[512];
            feff::common::logger().wlog("   Fixing edge energy from Elam table...");
            std::snprintf(slog, sizeof(slog), "   emu = %10.3f eV", emu * hart);
            feff::common::logger().wlog(slog);
        }
    }

    // Check logic for TDLDA/PMBSE flags
    if (ipmbse <= 0) itdlda = 0;
    if (nohole < 0) {
        if (ifxc != 0) {
            if (verbose) feff::common::logger().wlog(" Reset ifxc=0 since NOHOLE card is absent");
            ifxc = 0;
            if (ipmbse > 0) nonlocal = 0;
        }
        if (ipmbse == 3 && izstd == 0) {
            if (verbose) feff::common::logger().wlog(" Reset ipmbse=1 since NOHOLE card is absent");
            ipmbse = 1;
        }
    }
    if (izstd > 0 && itdlda > 0) {
        if (verbose) feff::common::logger().wlog(" Ignored PMBSE cards since TDLDA is present");
        itdlda = 0;
    }
    if (ipmbse == 2 && nonlocal > 0 && ifxc > 0) {
        if (verbose) feff::common::logger().wlog(" Reset ifxc=0 since core-hole potential is used.");
        ifxc = 0;
    }
    if (ipmbse == 1 && nonlocal > 0) nonlocal = 0;

    // Energy mesh
    double edge = xmu - vr0;
    if (!lopt) emu = emu - vr0;

    int ne, ne1, ne3;
    FeffComplex em[nex];
    double corr = 1.0;

    if (itdlda == 0) {
        phmesh(ipr2, ispec, edge, emu, vi0, gamach,
               xkmax, xkstep, vixan, ne, ne1, em, ik0, ne3);
    }

    // Atom r grid
    double ri[nrptx];
    for (int i = 0; i < 251; i++) {
        ri[i] = std::exp(-x0 + dx * i);
    }

    // Spin channels
    int nsp = 1;
    if (std::abs(ispin) == 1) nsp = nspx;

    // Scale spin density
    // (spinph is const, so we work with a local copy of dmag already)
    for (int iph = 0; iph <= nph; iph++) {
        for (int i = 0; i < 251; i++) {
            dmag[i + 251 * iph] *= spinph[iph];
        }
    }

    // Phase shift and cross-section arrays
    // Heap-allocated to avoid stack overflow (these total ~2 MB)
    FeffComplex eref[nex * nspx];
    std::vector<FeffComplex> ph(nex * (2 * ltot + 1) * nspx * (nphx + 1), FeffComplex(0.0, 0.0));
    int lmax_arr[nphx + 1];
    FeffComplex xsec_arr[nex * nspx];
    double xsnorm[nex * nspx];
    std::vector<FeffComplex> rkk(nex * 8 * nspx, FeffComplex(0.0, 0.0));

    // Working arrays — large ones on the heap to avoid stack overflow
    double vtotph[nrptx], rhoph[nrptx], dmagx[nrptx];
    double vvalph[nrptx], rhphvl[nrptx];
    double dgcx[nrptx], dpcx[nrptx];
    std::vector<double> dgcn_vec(nrptx * 30, 0.0);
    std::vector<double> dpcn_vec(nrptx * 30, 0.0);
    // Alias for code that uses dgcn/dpcn as 2D arrays [nrptx][30]
    auto dgcn = reinterpret_cast<double(*)[30]>(dgcn_vec.data());
    auto dpcn = reinterpret_cast<double(*)[30]>(dpcn_vec.data());
    double vjump;

    double xmuvr = xmu - vr0;

    for (int isp = 0; isp < nsp; isp++) {
        int ispinp;
        if (ispin != 0) {
            int idmag = (isp % 2 == 0) ? -1 : 1;
            if (nsp == 1) {
                idmag = (ispin < 0) ? -1 : 1;
            }
            // istprm would be called here for spin-dependent potential
            ispinp = ispin;
            if (std::abs(ispin) == 1 && nsp == 2) {
                ispinp = std::abs(ispin);
                if (isp == 0) ispinp = -ispinp;
            }
        } else {
            ispinp = 0;
        }

        // Calculate operators (szlz) if print level high enough
        if (ipr2 >= 3) {
            szlz(verbose, ispinp, ecv, nph, nat, rgrd, nohole, rfms2,
                 lfms2, lmaxph, edens, edenvl, dmag, vtot, vvalgs,
                 rmt, rnrm, ixc, rhoint, vint, xmuvr, jumprm,
                 xnval, iorb, x0, dx, xion, iunf, iz,
                 adgc, adpc, dgc, dpc, ihole,
                 rat, iphat, corr);
        }

        // Cross-section calculation for absorbing atom (iph=0)
        if (verbose) {
            feff::common::logger().wlog("    absorption cross section");
        }
        feff::common::fixvar(rmt[0], &edens[0], &vtot[0], &dmag[0],
                            vint, rhoint, dx, dxnew, jumprm,
                            vjump, ri, vtotph, rhoph, dmagx);
        feff::common::fixdsx(0, dx, dxnew, dgc, dpc,
                            reinterpret_cast<double*>(dgcn),
                            reinterpret_cast<double*>(dpcn));
        if (ixc % 10 >= 5) {
            int jrm2 = (jumprm > 0) ? 2 : jumprm;
            feff::common::fixvar(rmt[0], &edenvl[0], &vvalgs[0], &dmag[0],
                                vint, rhoint, dx, dxnew, jrm2,
                                vjump, ri, vvalph, rhphvl, dmagx);
            if (jumprm > 0) jumprm = 1;
        }
        int jnew;
        feff::common::fixdsp(dx, dxnew, dgc0_in, dpc0_in, dgcx, dpcx, jnew);

        if (itdlda == 0) {
            xsect(ipr2, dxnew, x0, ri, ne, ne1, ik0, em, edge,
                  ihole, emu, corr, dgcx, dpcx, jnew,
                  ixc0, lreal, rmt[0], rnrm[0], xmuvr, iPl,
                  vtotph, vvalph, rhoph, dmagx, rhphvl,
                  dgcn, dpcn,
                  reinterpret_cast<const double(*)[30]>(adgc),
                  reinterpret_cast<const double(*)[30]>(adpc),
                  &xsec_arr[isp * nex], &xsnorm[isp * nex],
                  reinterpret_cast<FeffComplex(*)[8]>(&rkk[isp * nex * 8]),
                  iz[0], xion[0], iunf, &xnval[0], izstd, &iorb[0], l2lp,
                  ipol, ispinp, le2, angks, ptz);
        }

        // Phase shifts for all unique potentials
        for (int iph = 0; iph <= nph; iph++) {
            if (verbose) {
                char slog[512];
                std::snprintf(slog, sizeof(slog),
                             "    phase shifts for unique potential%5d", iph);
                feff::common::logger().wlog(slog);
            }
            feff::common::fixvar(rmt[iph], &edens[251 * iph], &vtot[251 * iph],
                                &dmag[251 * iph], vint, rhoint, dx, dxnew, jumprm,
                                vjump, ri, vtotph, rhoph, dmagx);
            if (ixc % 10 >= 5) {
                int jrm2 = (jumprm > 0) ? 2 : jumprm;
                feff::common::fixvar(rmt[iph], &edenvl[251 * iph], &vvalgs[251 * iph],
                                    &dmag[251 * iph], vint, rhoint, dx, dxnew, jrm2,
                                    vjump, ri, vvalph, rhphvl, dmagx);
                if (jumprm > 0) jumprm = 1;
                feff::common::fixdsx(iph, dx, dxnew, dgc, dpc,
                                    reinterpret_cast<double*>(dgcn),
                                    reinterpret_cast<double*>(dpcn));
            }

            int itmp = (iph == 0) ? ihole : 0;

            // ph for this potential: ph[ie][-ltot:ltot] at offset [isp][iph]
            FeffComplex* ph_ptr = &ph[nex * (2 * ltot + 1) * (isp + nspx * iph)];

            // Fortran passes adgc(1,1,iph) — the subarray for potential iph.
            // Each iph slice occupies 10*30 = 300 doubles in the flat array.
            const double* adgc_iph = adgc + static_cast<size_t>(iph) * 10 * 30;
            const double* adpc_iph = adpc + static_cast<size_t>(iph) * 10 * 30;

            phase(iph, dxnew, x0, ri, ne, ne1, ne3, em, ixc, nsp,
                  lreal, rmt[iph], rnrm[iph], xmuvr, iPl,
                  vtotph, vvalph, rhoph, dmagx, rhphvl,
                  dgcn, dpcn,
                  reinterpret_cast<const double(*)[30]>(adgc_iph),
                  reinterpret_cast<const double(*)[30]>(adpc_iph),
                  &eref[isp * nex], ph_ptr, lmax_arr[iph],
                  iz[iph], itmp, xion[iph], iunf, &xnval[30 * iph], ispinp);
        }
    } // spin loop

    // Prepare output columns
    double col1[nex], col2[nex], col3[nex], col4[nex], col5[nex];
    if (std::abs(ispin) != 1 || nspx == 1) {
        for (int ie = 0; ie < ne; ie++) {
            col1[ie] = std::real(em[ie]) * hart;
            col2[ie] = std::imag(em[ie]) * hart;
            col3[ie] = xsnorm[ie];
            col4[ie] = std::real(xsec_arr[ie]);
            col5[ie] = std::imag(xsec_arr[ie]);
        }
    } else {
        // Spin-polarized case: average
        for (int ie = 0; ie < ne; ie++) {
            col1[ie] = std::real(em[ie]) * hart;
            col2[ie] = std::imag(em[ie]) * hart;
            col3[ie] = (xsnorm[ie] + xsnorm[ie + nex * (nspx - 1)]) / 2.0;
            col4[ie] = std::real(xsec_arr[ie] + xsec_arr[ie + nex * (nspx - 1)]);
            col5[ie] = std::imag(xsec_arr[ie] + xsec_arr[ie + nex * (nspx - 1)]);

            // Normalize rkk
            // rkk is stored in C row-major: rkk[isp * nex * 8 + ie * 8 + k]
            double sum = xsnorm[ie] + xsnorm[ie + nex * (nspx - 1)];
            if (sum > 0.0) {
                double xnorm1 = std::sqrt(2.0 * xsnorm[ie] / sum);
                double xnorm2 = std::sqrt(2.0 * xsnorm[ie + nex * (nspx - 1)] / sum);
                for (int k = 0; k < 8; k++) {
                    rkk[ie * 8 + k] *= xnorm1;
                    rkk[(nspx - 1) * nex * 8 + ie * 8 + k] *= xnorm2;
                }
            }
        }
    }

    // NaN check
    for (int ie = 0; ie < ne; ie++) {
        if (std::isnan(col3[ie])) col3[ie] = 0.0;
        if (std::isnan(col4[ie])) col4[ie] = 0.0;
        if (std::isnan(col5[ie])) col5[ie] = 0.0;
    }

    // Write xsect.json (cross-section data for FF2X)
    if (wrxsec) {
        // Convert char[][80] titles to std::string[]
        std::vector<std::string> stitle(ntitle);
        for (int i = 0; i < ntitle; ++i) {
            stitle[i] = std::string(title[i], 80);
            // Trim trailing spaces
            size_t end = stitle[i].find_last_not_of(' ');
            if (end != std::string::npos) stitle[i].resize(end + 1);
            else stitle[i].clear();
        }
        feff::json_io::write_xsect_json(ntitle, stitle.data(), s02, erelax,
                                         wp, edge, emu,
                                         gamach * hart, ne, ne1, ik0,
                                         col1, col2, col3, col4, col5);
        if (verbose) {
            feff::common::logger().wlog("    Wrote xsect.json");
        }
    }

    // Write diagnostic phase files
    if (ipr2 >= 2) {
        wphase(nph, em, eref, lmax_arr, ne, ph.data(), ntitle, title);
    }

    // Write phase.pad
    wrxsph(phpad, nsp, ne, ne1, ne3, nph, ihole, rnrmav, xmuvr,
           edge, ik0, ixc, rs, vint, em, eref, lmax_arr, iz, potlbl, ph.data(), rkk.data());

    // Calculate AXAFS
    if (ipr2 >= 1) {
        axafs(em, emu, xsec_arr, ne1, ik0);
    }
}

} // namespace feff::xsph
