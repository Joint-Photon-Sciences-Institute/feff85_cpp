// Find MT radii and interstitial parameters.
// Converted from src/POT/istprm.f

#include "istprm.hpp"
#include "sidx.hpp"
#include "fermi.hpp"
#include "../exch/vbh.hpp"
#include "../exch/edp.hpp"
#include "../math/distance.hpp"
#include "../common/radial_grid.hpp"
#include "../common/logging.hpp"
#include <feff/constants.hpp>
#include <cmath>
#include <sstream>
#include <complex>

// Include actual headers for movrlp and ovp2mt
#include "movrlp.hpp"
#include "ovp2mt.hpp"

namespace feff::pot {

double calcvl(double r1, double r2, double r)
{
    double xl = (r1 * r1 - r2 * r2 + r * r) / (2.0 * r);
    double h = r1 - xl;
    return (pi / 3.0) * h * h * (3.0 * r1 - h);
}

void istprm(int nph, int nat, const int* iphat, const double* rat,
            const int* iatph, const double* xnatph,
            const int* novr, const int* iphovr, const int* nnovr, const double* rovr,
            const double* folp, double* folpx, int iafolp,
            double* edens, double* edenvl, int idmag,
            double* dmag, double* vclap, double* vtot, double* vvalgs,
            int* imt, int* inrm, double* rmt, double* rnrm,
            int ixc, double& rhoint, double& vint, double& rs, double& xf,
            double xmu, double& xmunew,
            double& rnrmav, double& qtotel, int& inters, double totvol)
{
    using feff::common::ii;
    using feff::common::rr;
    using feff::common::xx;

    constexpr int s251 = 251;
    constexpr double big = 5000.0;
    constexpr int novp = 40;

    // Work space
    double ri[251];
    static bool lnear[nphx + 1];
    static int inn[nphx + 1];
    static double rnnmin[nphx + 1];

    // Overlap matrix workspace
    constexpr int cmovp_dim = novp * (nphx + 1) + 1;
    std::complex<float> cmovp[cmovp_dim * cmovp_dim];
    int ipiv[cmovp_dim];

    // Find rmt from rnrm only on first call of istprm (rmt[0] <= 0)
    if (rmt[0] <= 0.0) {
        for (int iph = 0; iph <= nph; ++iph) {
            lnear[iph] = false;
        }

        for (int iph = 0; iph <= nph; ++iph) {
            double voltot = 0.0;
            double rmtavg = 0.0;
            inrm[iph] = ii(rnrm[iph]);

            if (novr[iph] > 0) {
                // Overlap explicitly defined by overlap card
                double rnear = big;
                inters = inters % 6;

                for (int iovr = 0; iovr < novr[iph]; ++iovr) {
                    double rnn = rovr[iph * novrx + iovr];
                    int inph = iphovr[iph * novrx + iovr];
                    if (rnn <= rnear) {
                        rnear = rnn;
                        rnnmin[iph] = rnn;
                        inn[iph] = inph;
                    }
                    // Don't avg if Norman spheres don't overlap
                    if (rnrm[iph] + rnrm[inph] > rnn) {
                        double voltmp = calcvl(rnrm[iph], rnrm[inph], rnn);
                        voltmp = voltmp + calcvl(rnrm[inph], rnrm[iph], rnn);
                        double rmttmp = rnn * folp[iph] * rnrm[iph] /
                                        (rnrm[iph] + rnrm[inph]);
                        int ntmp = nnovr[iph * novrx + iovr];
                        rmtavg = rmtavg + rmttmp * voltmp * ntmp;
                        voltot = voltot + voltmp * ntmp;
                    }
                }
            } else {
                int iat = iatph[iph];  // 1-based
                double rnear = big;
                rmt[iph] = big;

                for (int inat = 0; inat < nat; ++inat) {
                    if (inat == iat - 1) continue;  // iat is 1-based
                    double rnn = feff::math::dist(&rat[inat * 3], &rat[(iat - 1) * 3]);
                    int inph = iphat[inat];
                    if (rnn <= rnear) {
                        rnear = rnn;
                        rnnmin[iph] = rnn;
                        inn[iph] = inph;
                    }
                    // Don't avg if Norman spheres don't overlap
                    if (rnrm[iph] + rnrm[inph] >= rnn) {
                        if (inters < 6) {
                            // Norman prescription
                            double voltmp = calcvl(rnrm[iph], rnrm[inph], rnn);
                            voltmp = voltmp + calcvl(rnrm[inph], rnrm[iph], rnn);
                            double rmttmp = rnn * folp[iph] * rnrm[iph] /
                                            (rnrm[iph] + rnrm[inph]);
                            rmtavg = rmtavg + rmttmp * voltmp;
                            voltot = voltot + voltmp;
                        } else {
                            // Matching point prescription
                            for (int i = inrm[iph]; i >= 1; --i) {
                                int j = ii(rnn - rnrm[iph]);
                                if (vclap[iph * s251 + (i - 1)] <=
                                    vclap[inph * s251 + (j - 1)]) {
                                    double d1 = (vclap[iph * s251 + i] -
                                                 vclap[iph * s251 + (i - 1)]) /
                                                (rr(i + 1) - rr(i));
                                    double d2 = (vclap[inph * s251 + (j - 1)] -
                                                 vclap[inph * s251 + (j - 2)]) /
                                                (rr(j) - rr(j - 1));
                                    rmtavg = rr(i) +
                                             (vclap[inph * s251 + (j - 1)] +
                                              d2 * (rnn - rr(i) - rr(j)) -
                                              vclap[iph * s251 + (i - 1)]) /
                                             (d1 + d2);
                                    break;
                                }
                            }
                            if (rmtavg < rmt[iph]) rmt[iph] = rmtavg;
                        }
                    }
                }
            }

            // Special situation if rnrm is too close or larger than NN distance
            if (novr[iph] > 0 || inters < 6) {
                double rnear = rnnmin[iph];
                if (rnrm[iph] >= rnear) lnear[iph] = true;
            }

            if (rmtavg <= 0.0) {
                int iat = iatph[iph];
                std::ostringstream slog;
                slog << " WARNING: NO ATOMS CLOSE ENOUGH TO OVERLAP ATOM"
                     << iat << ",  UNIQUE POT" << iph
                     << "!!  Rmt set to Rnorman.  May be error in input file.";
                feff::common::logger().wlog(slog.str());
                rmt[iph] = rnrm[iph];
            } else if (inters < 6) {
                // Norman prescription
                rmt[iph] = rmtavg / voltot;
                double rnear = rnnmin[iph];
                if (rmt[iph] >= rnear) {
                    feff::common::logger().wlog(
                        " Rmt >= distance to nearest neighbor.  Not physically meaningful.");
                    feff::common::logger().wlog(
                        " FEFF may crash.  Look for error in ATOM list or OVERLAP cards.");
                }
                if (rnrm[iph] >= rnear) {
                    int imax = ii(rnear) - 1;
                    // Until loop: find where potential decreases
                    while (vclap[iph * s251 + (imax - 1)] <
                           vclap[iph * s251 + imax]) {
                        imax = imax - 1;
                    }
                    rmt[iph] = std::exp(xx(imax)) - 0.0001;
                }
            }
        }

        // Set maximum value for folp(iph) if AFOLP is in use
        for (int iph = 0; iph <= nph; ++iph) {
            double temp;
            if (iafolp > 0) {
                temp = 0.2 + 0.8 * rnrm[iph] / rmt[iph];
            } else {
                temp = 0.3 + 0.7 * rnrm[iph] / rmt[iph];
            }
            if (temp < folpx[iph]) folpx[iph] = temp;
            temp = rnnmin[iph] / rmt[iph] / 1.06;
            if (temp < folpx[iph]) folpx[iph] = temp;
            temp = std::exp(-(novp - 3) * 0.05);

            if (lnear[iph]) {
                temp = rnnmin[iph] / (rmt[iph] * 1.05 + temp * rmt[inn[iph]]);
                if (temp < folpx[iph]) folpx[iph] = temp;
                if (temp < folpx[inn[iph]]) folpx[inn[iph]] = temp;
            } else {
                temp = (rnnmin[iph] - rnrm[iph]) / (temp * rmt[inn[iph]]);
                if (temp < folpx[inn[iph]]) folpx[inn[iph]] = temp;
            }
        }
    }
    // End of finding rmt from rnrm on first call of istprm.

    // Need potential with ground state xc, put it into vtot
    for (int iph = 0; iph <= nph; ++iph) {
        int imax_local;
        sidx(&edens[iph * s251], 250, rmt[iph], rnrm[iph],
             imax_local, imt[iph], inrm[iph]);

        for (int i = 0; i < imax_local; ++i) {
            // i is 0-based; Fortran i is 1-based
            double rs_local;
            double xmag;
            if (edens[iph * s251 + i] <= 0.0) {
                rs_local = 100.0;
                xmag = 1.0;
            } else {
                rs_local = std::pow(edens[iph * s251 + i] / 3.0, -third);
                xmag = 1.0 + idmag * dmag[iph * s251 + i];
            }

            double vvbh;
            feff::exch::vbh(rs_local, xmag, vvbh);
            vtot[iph * s251 + i] = vclap[iph * s251 + i] + vvbh;

            if ((ixc % 10) == 5) {
                double rsval = 10.0;
                if (edenvl[iph * s251 + i] > 0.00001)
                    rsval = std::pow(edenvl[iph * s251 + i] / 3.0, -third);
                if (rsval > 10.0) rsval = 10.0;
                double xmagvl = 1.0 + idmag * dmag[iph * s251 + i] *
                                edens[iph * s251 + i] / edenvl[iph * s251 + i];
                double vvbhvl;
                feff::exch::vbh(rsval, xmagvl, vvbhvl);
                vvalgs[iph * s251 + i] = vclap[iph * s251 + i] + vvbhvl;
            } else if ((ixc % 10) >= 6) {
                double rscore;
                if (edens[iph * s251 + i] <= edenvl[iph * s251 + i]) {
                    rscore = 101.0;
                } else {
                    rscore = std::pow((edens[iph * s251 + i] - edenvl[iph * s251 + i]) / 3.0,
                                      -third);
                }
                double rsmag = std::pow(edens[iph * s251 + i] *
                               (1.0 + idmag * dmag[iph * s251 + i]) / 3.0, -third);
                double xfmag = fa / rsmag;
                double vrdh;
                feff::exch::edp(rscore, xfmag, vrdh);
                vvalgs[iph * s251 + i] = vclap[iph * s251 + i] + vvbh - vrdh;
            } else {
                vvalgs[iph * s251 + i] = 0.0;
            }
        }
    }

    // Calculate average Norman radius and interstitial volume
    rnrmav = 0.0;
    double xn = 0.0;
    double volint = 0.0;
    for (int iph = 0; iph <= nph; ++iph) {
        rnrmav = rnrmav + xnatph[iph] * rnrm[iph] * rnrm[iph] * rnrm[iph];
        volint = volint - xnatph[iph] * rmt[iph] * rmt[iph] * rmt[iph];
        xn = xn + xnatph[iph];
    }
    if (totvol <= 0.0) {
        volint = 4.0 * pi / 3.0 * (volint + rnrmav);
    } else {
        volint = 4.0 * pi / 3.0 * volint + totvol;
    }
    rnrmav = std::pow(rnrmav / xn, third);

    rs = 0.0;
    vint = 0.0;
    rhoint = 0.0;

    movrlp(nph, nat, iphat, rat, iatph, xnatph,
           novr, iphovr, nnovr, rovr,
           imt, rmt, rnrm, ri, lnear,
           cmovp, ipiv, volint, inters);

    // If no contribution to interstitial from any atom, die
    if (volint <= 0.0) {
        feff::common::logger().wlog(" No interstitial density.  Check input file.");
        throw std::runtime_error("ISTPRM");
    }

    // Find interstitial density
    ovp2mt(nph, edens, 0, qtotel, ri, xnatph, lnear,
           inrm, imt, rnrm, rmt, cmovp, ipiv, rhoint, inters);
    rhoint = 4.0 * pi * rhoint / volint;

    if (ixc >= 5) {
        double vintvl;
        ovp2mt(nph, vvalgs, 1, qtotel, ri, xnatph, lnear,
               inrm, imt, rnrm, rmt, cmovp, ipiv, vintvl, inters);
    }

    // Find potential inside MT sphere and vint
    ovp2mt(nph, vtot, 1, qtotel, ri, xnatph, lnear,
           inrm, imt, rnrm, rmt, cmovp, ipiv, vint, inters);

    if (vint >= xmu) {
        feff::common::logger().wlog(
            " WARNING:interstitial level found above Fermi level");
        feff::common::logger().wlog(
            "  Results may be unreliable. See manual for details");
        vint = xmu - 0.05;
        ovp2mt(nph, vtot, 2, qtotel, ri, xnatph, lnear,
               inrm, imt, rnrm, rmt, cmovp, ipiv, vint, inters);
    }

    fermi(rhoint, vint, xmunew, rs, xf);
}

} // namespace feff::pot
