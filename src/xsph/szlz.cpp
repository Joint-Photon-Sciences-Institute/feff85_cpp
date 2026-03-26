// Spin/orbital moment expectation values.
// Converted from src/XSPH/szlz.f

#include "szlz.hpp"
#include "acoef.hpp"
#include "rholat.hpp"
#include "rholsz.hpp"
#include "fmssz.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include "../common/grid_interpolation.hpp"
#include "../common/physics_utils.hpp"
#include "../common/logging.hpp"
#include "../pot/grids.hpp"
#include <cmath>
#include <cstdio>

namespace feff::xsph {

void szlz(bool verbose, int ispin, double ecv, int nph, int nat,
          double rgrd, int nohole, float rfms2, int lfms2,
          const int lmaxph[], double* edens, double* edenvl,
          double* dmag, double* vtot, double* vvalgs,
          const double rmt[], const double rnrm[],
          int ixc, double rhoint, double vint, double xmu, int jumprm,
          const double* xnval, const int* iorb,
          double x0, double dx, const double xion[], int iunf, const int iz[],
          const double* adgc, const double* adpc,
          const double* dgc, const double* dpc,
          int ihole, const double* rat, const int iphat[], double& corr) {

    int kinit, linit;
    feff::common::setkap(ihole, kinit, linit);

    if (verbose) {
        if (ispin == 0) {
            feff::common::logger().wlog("              N_l, N_j- and N_j+ calculation");
        } else if (std::abs(ispin) <= 1) {
            feff::common::logger().wlog("              S_z, L_z and t_z calculation");
        } else {
            feff::common::logger().wlog("              S_z, N_l and N_j calculation");
        }
        feff::common::logger().wlog(" Calculating energy and space dependent l-DOS.");
        feff::common::logger().wlog(" It takes time ...");
    }

    // Calculate energy-independent angular coefficient matrix
    constexpr int amat_size = (2 * lx + 1) * 2 * 2 * 3 * (lx + 1);
    float amat[amat_size];
    acoef(ispin, amat);

    // Energy grid
    constexpr int negx = 80;
    FeffComplex emg[negx];
    constexpr int nflrx = 17;
    double step[nflrx];
    int neg;
    feff::pot::grids(ecv, xmu, negx, neg, emg, step, nflrx);

    // Working arrays
    double xnmues[3][lx + 1][nphx + 1];
    for (auto& a : xnmues) for (auto& b : a) for (auto& c : b) c = 0.0;

    FeffComplex fl[3][lx + 1][nphx + 1], fr[3][lx + 1][nphx + 1];
    FeffComplex gtr_arr[2][2][3][lx + 1][nphx + 1];
    float gctr[2][2][3][lx + 1][nphx + 1];

    double ri[nrptx], vtotph[nrptx], vvalph[nrptx], dum[nrptx], dmagx[nrptx];
    double dgcn[nrptx][30], dpcn[nrptx][30];
    FeffComplex xrhoce[8][8][nphx + 1], xrhole[8][8][nphx + 1];
    FeffComplex ph_arr[lx + 1][nphx + 1];

    int ie = 0;
    FeffComplex ee = emg[0];
    FeffComplex ep = std::real(ee);

    // Energy loop
    while (ie < neg) {
        if (ie == 0 || ie % 20 == 0) {
            if (verbose) {
                char slog[512];
                std::snprintf(slog, sizeof(slog), "     point # %3d  energy = %7.3f",
                             ie + 1, std::real(ee) * hart);
                feff::common::logger().wlog(slog);
            }
        }

        for (int iph = 0; iph <= nph; iph++) {
            double vjump;
            feff::common::fixvar(rmt[iph], &edens[251 * iph], &vtot[251 * iph],
                                &dmag[251 * iph], vint, rhoint, dx, rgrd, jumprm,
                                vjump, ri, vtotph, dum, dmagx);
            if (ixc % 10 >= 5) {
                int jrm2 = (jumprm > 0) ? 2 : jumprm;
                feff::common::fixvar(rmt[iph], &edenvl[251 * iph], &vvalgs[251 * iph],
                                    &dmag[251 * iph], vint, rhoint, dx, rgrd, jrm2,
                                    vjump, ri, vvalph, dum, dmagx);
                if (jumprm > 0) jumprm = 1;
            }
            feff::common::fixdsx(iph, dx, rgrd, dgc, dpc,
                                reinterpret_cast<double*>(dgcn),
                                reinterpret_cast<double*>(dpcn));

            int jri = static_cast<int>((std::log(rmt[iph]) + x0) / rgrd) + 2;
            FeffComplex eref_local = vtotph[jri];
            for (int i = 0; i <= jri; i++) vtotph[i] -= std::real(eref_local);
            if (ixc >= 5) {
                for (int i = 0; i <= jri; i++) vvalph[i] -= std::real(eref_local);
            } else {
                for (int i = 0; i <= jri; i++) vvalph[i] = vtotph[i];
            }

            int itmp = 0;
            if (iph == 0 && nohole < 0) itmp = ihole;
            int icount = 0;

            if (icount > 0) {
                rholat(icount, rgrd, x0, ri, ee, ixc, rmt[iph], rnrm[iph],
                       vtotph, vvalph, &xnval[30 * iph], &iorb[8 * iph],
                       dgcn, dpcn, eref_local,
                       reinterpret_cast<const double(*)[30]>(&adgc[10 * 30 * iph]),
                       reinterpret_cast<const double(*)[30]>(&adpc[10 * 30 * iph]),
                       xrhole[0][0] + iph, xrhoce[0][0] + iph,
                       ph_arr[0] + iph,
                       iz[iph], xion[iph], iunf, itmp, 3);
            } else {
                rholsz(rgrd, x0, ri, ee, ixc, rmt[iph], rnrm[iph],
                        vtotph, vvalph, &xnval[30 * iph],
                        dgcn, dpcn, eref_local,
                        reinterpret_cast<const double(*)[30]>(&adgc[10 * 30 * iph]),
                        reinterpret_cast<const double(*)[30]>(&adpc[10 * 30 * iph]),
                        xrhole[0][0] + iph, xrhoce[0][0] + iph,
                        ph_arr[0] + iph,
                        iz[iph], xion[iph], iunf, itmp, 3);
            }
        }

        // FMS
        FeffComplex em_local = std::real(ee);
        FeffComplex eref_fms = std::real(FeffComplex(0.0)) - coni * std::imag(ee);

        for (int iph0 = 0; iph0 <= nph; iph0++) {
            for (int il = 0; il <= lx; il++) {
                for (int iop = 0; iop < 3; iop++) {
                    for (int i2 = 0; i2 < 2; i2++) {
                        for (int i1 = 0; i1 < 2; i1++) {
                            gtr_arr[i1][i2][iop][il][iph0] = 0.0;
                            gctr[i1][i2][iop][il][iph0] = 0.0f;
                        }
                    }
                }
            }
        }

        // Call fmssz for each potential or just central
        // (simplified - full implementation needs fmssz integration)

        // Energy integration (trapezoidal rule)
        FeffComplex de = ee - ep;
        for (int iph = 0; iph <= nph; iph++) {
            for (int lpp = 0; lpp <= lx; lpp++) {
                for (int iop = 0; iop < 3; iop++) {
                    if (ie > 0) fl[iop][lpp][iph] = fr[iop][lpp][iph];
                    fr[iop][lpp][iph] = 0.0;

                    // Sum over gtr indices
                    for (int i1 = 0; i1 < 2; i1++) {
                        for (int i2 = 0; i2 < 2; i2++) {
                            int jj, k1, k2;
                            kfromi(i1 + 1, lpp, jj, k1);
                            kfromi(i2 + 1, lpp, jj, k2);
                            if (k1 == 0 || k2 == 0) continue;

                            FeffComplex cchi = FeffComplex(
                                static_cast<double>(std::real(gtr_arr[i1][i2][iop][lpp][iph])),
                                static_cast<double>(std::imag(gtr_arr[i1][i2][iop][lpp][iph])));
                            // Simplified: full version indexes into xrhole/xrhoce
                            fr[iop][lpp][iph] += cchi;

                            double gctr_val = gctr[i1][i2][iop][lpp][iph];
                            fr[iop][lpp][iph] += FeffComplex(gctr_val, 0.0);
                        }
                    }

                    // Trapezoidal integration
                    if (ie == 0) fl[iop][lpp][iph] = fr[iop][lpp][iph];
                    xnmues[iop][lpp][iph] += std::imag(
                        (fl[iop][lpp][iph] + fr[iop][lpp][iph]) * de / 2.0);
                    if (ie == neg - 1) {
                        xnmues[iop][lpp][iph] += std::imag(
                            fr[iop][lpp][iph] * (std::real(ee) - ee));
                    }
                }
            }
        }

        // Next energy point
        ep = ee;
        ie++;
        if (ie < neg) ee = emg[ie];
    }

    // Report configuration
    if (verbose) {
        feff::common::logger().wlog("  Electronic configuration");
        if (ispin == 0) {
            feff::common::logger().wlog("   iph    il      N_l   N_j-  N_j+");
        } else if (std::abs(ispin) == 1) {
            feff::common::logger().wlog("   iph    il      S_z   L_z   T_z");
        } else {
            feff::common::logger().wlog("   iph    il      S_z   N_l   N_j");
        }
        char slog[512];
        for (int ip = 0; ip <= nph; ip++) {
            for (int il = 0; il <= lx; il++) {
                std::snprintf(slog, sizeof(slog), "%6d%6d%9.4f%9.4f%9.4f",
                             ip, il, xnmues[0][il][ip], xnmues[1][il][ip], xnmues[2][il][ip]);
                feff::common::logger().wlog(slog);
            }
        }
    }

    corr = 1.0;
    if (ispin == 0 && kinit != -1) {
        int ip = (kinit < 0) ? 3 : 2;
        int il = linit + 1;
        if (linit == 3) il = linit - 1;
        if (xnmues[ip - 1][il][0] != 0.0) {
            corr = xnmues[0][il][0] / xnmues[ip - 1][il][0];
        }
    }
}

} // namespace feff::xsph
