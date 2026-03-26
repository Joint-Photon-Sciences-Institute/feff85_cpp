// Overlap matrix construction for muffin-tin potential decomposition.
// Converted from: src/POT/movrlp.f
//
// Constructs the overlap matrix based on geometry of overlapped
// muffin-tin spheres, then performs LU decomposition for use by ovp2mt.

#include "movrlp.hpp"

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include "../common/logging.hpp"
#include "../par/parallel.hpp"
#include "../math/distance.hpp"

#include <cmath>
#include <complex>
#include <cstring>

// Forward declarations for LAPACK
extern "C" {
    // LAPACK: LU factorization (single-precision complex)
    void cgetrf_(const int* m, const int* n, std::complex<float>* a,
                 const int* lda, int* ipiv, int* info);
}

namespace feff::pot {

// Forward declaration: volume of intersection of two spheres
// calcvl(r1, r2, rnn) - volume cap of sphere r1 cut by sphere r2 at distance rnn
double calcvl(double r1, double r2, double rnn);

// Helper: convert radius to grid index (Fortran function ii)
// ii(r) = int((log(r) + 8.8) / 0.05) + 1  (1-based Fortran index)
static int ii_func(double r) {
    if (r <= 0.0) return 1;
    return static_cast<int>((std::log(r) + 8.8) / 0.05) + 1;
}

void movrlp(int nph, int nat, const int* iphat, const double* rat,
            const int* iatph, const double* xnatph,
            const int* novr, const int* iphovr, const int* nnovr,
            const double* rovr,
            const int* imt, const double* rmt, const double* rnrm,
            double* ri, bool* lnear,
            std::complex<float>* cmovp, int* ipiv,
            double& volint, int inters)
{
    using CFloat = std::complex<float>;
    auto& log = feff::common::logger();
    char slog[512];

    int ntmp = 0;
    int iat0 = -999;

    // Extract ipot and irav from inters
    int ipot = inters % 2;
    int irav = (inters - ipot) / 2;

    // Initialize ri (Loucks grid)
    for (int i = 0; i < 251; i++) {
        ri[i] = std::exp(-8.8 + i * 0.05);
    }
    double exphx = std::exp(0.025);

    // Initialize cmovp as unit matrix up to ncp
    int ncp = novp * (nph + 1) + 1;

    // cmovp is stored as column-major [istatx_dim * istatx_dim]
    // Access: cmovp[i1 + istatx_dim * i2]  (0-based)
    for (int i2 = 0; i2 < ncp; i2++) {
        for (int i1 = 0; i1 < ncp; i1++) {
            cmovp[i1 + istatx_dim * i2] = CFloat(0.0f, 0.0f);
            if (i1 == i2) cmovp[i1 + istatx_dim * i2] = CFloat(1.0f, 0.0f);
            if (i2 == ncp - 1) cmovp[i1 + istatx_dim * i2] = CFloat(0.01f, 0.0f);
        }
    }

    // bmat(nphx+1, novp*(nphx+1)) - stored as bmat[ix1 + (nphx+1)*ix2]
    float bmat[(nphx + 1) * novp * (nphx + 1)];
    std::memset(bmat, 0, sizeof(bmat));

    double xn = 0.0;

    // Main loop over unique potentials (label 200)
    for (int ip1 = 0; ip1 <= nph; ip1++) {
        int nlast;
        if (novr[ip1] > 0) {
            nlast = novr[ip1];
        } else {
            iat0 = iatph[ip1];
            ntmp = 1;
            nlast = nat;
        }

        double rav;
        if (irav == 1) {
            rav = (rmt[ip1] + rnrm[ip1]) / 2.0;
        } else if (irav == 0) {
            rav = rnrm[ip1];
        } else {
            // ri[imt[ip1]] corresponds to Fortran ri(imt(ip1)+1) since imt is 1-based
            rav = ri[imt[ip1]];
        }
        if (lnear[ip1]) rav = ri[imt[ip1]];

        // Loop over neighbors (label 190)
        for (int iat = 1; iat <= nlast; iat++) {
            int ip2;
            double rnn;
            if (novr[ip1] > 0) {
                // Fortran: nnovr(iat, ip1), iphovr(iat, ip1), rovr(iat, ip1)
                // These are (novrx, 0:nphx) arrays
                ntmp = nnovr[(iat - 1) + novrx * ip1];
                ip2  = iphovr[(iat - 1) + novrx * ip1];
                rnn  = rovr[(iat - 1) + novrx * ip1];
            } else {
                if (iat == iat0) continue; // goto 190
                ip2 = iphat[iat - 1];  // 0-based atom index
                // rat is (3, natx) column-major
                // dist(rat(1,iat0), rat(1,iat)) -> distance between atoms iat0 and iat
                rnn = feff::math::dist(&rat[3 * (iat0 - 1)], &rat[3 * (iat - 1)]);
            }

            // Correct for double counting volume and area
            if (rnn < rmt[ip1] + rmt[ip2]) {
                // Note: Fortran original calls calcvl(rmt(ip1), rmt(ip2), rnn) twice
                // (not with swapped args). This is a Fortran bug, but we reproduce it
                // to match Fortran results. See movrlp.f lines 96-98.
                volint += xnatph[ip1] * ntmp *
                    (calcvl(rmt[ip1], rmt[ip2], rnn) +
                     calcvl(rmt[ip1], rmt[ip2], rnn)) / 2.0;
            }

            // Fill bmat: expression for vtot(jri)
            int ix1 = ip1;  // 0-based row in bmat (Fortran ip1+1, 1-based)

            if (rav + rmt[ip2] <= rnn) goto label_100;

            {
                int imin2_check = ii_func(rnn - rav);
                if (imt[ip2] - imin2_check >= novp - 1) {
                    std::snprintf(slog, sizeof(slog),
                        " FOLP for POTENTIAL type %3d is too big.", ip1);
                    log.wlog(slog);
                    log.wlog(" Reduce overlap using FOLP and rerun");
                    feff::par::par_stop("MOVRLP-1");
                }
            }

            {
                int imin2 = imt[ip2] - novp + 1;  // 1-based Fortran index

                // Loop over i2 from imin2 to imt(ip2) (1-based)
                for (int i2 = imin2; i2 <= imt[ip2]; i2++) {
                    int i2_0 = i2 - 1;  // 0-based grid index
                    double r1 = ri[i2_0] / exphx;
                    double r2 = ri[i2_0] * exphx;
                    if (i2 == imt[ip2]) r2 = rmt[ip2];
                    if (i2 == imt[ip2]) r1 = (r1 + 2.0 * ri[imt[ip2] - 1] - rmt[ip2]) / 2.0;
                    if (i2 == imt[ip2] - 1) r2 = (r2 + 2.0 * ri[imt[ip2] - 1] - rmt[ip2]) / 2.0;
                    if (r2 + rav < rnn) continue;  // goto 80

                    if (r1 + rav < rnn) {
                        // Linear interpolation
                        double xr = (rnn - rav - r1) / (r2 - r1);
                        r1 = rnn - rav;
                        double temp = (r2 * r2 - r1 * r1) / (4.0 * rnn * rav) * ntmp;
                        int ind2 = i2 + 1;
                        if (i2 == imt[ip2]) ind2 = i2 - 1;
                        int ind2_0 = ind2 - 1;
                        xr = xr * (r2 - ri[i2_0]) / (ri[ind2_0] - ri[i2_0]);

                        // bmat indices: bmat[ix1 + (nphx+1) * ix2]
                        int ix2 = ip2 * novp + (i2 - imin2);
                        bmat[ix1 + (nphx + 1) * ix2] += static_cast<float>(temp * (1.0 - xr));
                        ix2 = ip2 * novp + (ind2 - imin2);
                        bmat[ix1 + (nphx + 1) * ix2] += static_cast<float>(temp * xr);
                    } else {
                        double temp = (r2 * r2 - r1 * r1) / (4.0 * rnn * rav) * ntmp;
                        int ix2 = ip2 * novp + (i2 - imin2);
                        bmat[ix1 + (nphx + 1) * ix2] += static_cast<float>(temp);
                    }
                } // end i2 loop (label 80)
            }

label_100:
            // Fill cmovp: expression for vtot(i), i < jri
            if (rmt[ip1] + rmt[ip2] <= rnn) continue;  // goto 190

            {
                int imin1 = imt[ip1] - novp + 1;  // 1-based
                int imin2 = imt[ip2] - novp + 1;  // 1-based

                {
                    int imin1_check = ii_func(rnn - rmt[ip2]);
                    int imin2_check = ii_func(rnn - rmt[ip1]);
                    if (imt[ip1] - imin1_check >= novp - 1 ||
                        imt[ip2] - imin2_check >= novp - 1) {
                        feff::par::par_stop("tell authors to INCREASE NOVP");
                    }
                }

                // Loop over i1 from imin1 to imt(ip1) (label 180)
                for (int i1 = imin1; i1 <= imt[ip1]; i1++) {
                    int i1_0 = i1 - 1;
                    double ri1 = ri[i1_0] / exphx;
                    double ri2 = ri[i1_0] * exphx;
                    if (i1 == imt[ip1]) ri2 = rmt[ip1];
                    if (i1 == imt[ip1]) ri1 = (ri1 + 2.0 * ri[imt[ip1] - 1] - rmt[ip1]) / 2.0;
                    if (i1 == imt[ip1] - 1) ri2 = (ri2 + 2.0 * ri[imt[ip1] - 1] - rmt[ip1]) / 2.0;

                    // cmovp index for i1: ix1 = (i1 - imin1) + ip1*novp  (0-based)
                    int ix1_cm = (i1 - imin1) + ip1 * novp;

                    // Loop over i2 from imin2 to imt(ip2) (label 170)
                    for (int i2 = imin2; i2 <= imt[ip2]; i2++) {
                        int i2_0 = i2 - 1;
                        double r1 = ri[i2_0] / exphx;
                        double r2 = ri[i2_0] * exphx;
                        if (i2 == imt[ip2]) r2 = rmt[ip2];
                        if (i2 == imt[ip2]) r1 = (r1 + 2.0 * ri[imt[ip2] - 1] - rmt[ip2]) / 2.0;
                        if (i2 == imt[ip2] - 1) r2 = (r2 + 2.0 * ri[imt[ip2] - 1] - rmt[ip2]) / 2.0;

                        if (r2 + ri2 < rnn) continue;  // goto 170

                        // Calculate volume of intersection
                        double temp = calcvl(ri2, r2, rnn) + calcvl(r2, ri2, rnn);
                        if (ri1 + r2 > rnn)
                            temp -= calcvl(ri1, r2, rnn) + calcvl(r2, ri1, rnn);
                        if (ri2 + r1 > rnn)
                            temp -= calcvl(ri2, r1, rnn) + calcvl(r1, ri2, rnn);
                        if (ri1 + r1 > rnn)
                            temp += calcvl(ri1, r1, rnn) + calcvl(r1, ri1, rnn);

                        // Divide by shell volume
                        temp = temp / (4.0 / 3.0 * pi * (ri2 * ri2 * ri2 - ri1 * ri1 * ri1)) * ntmp;

                        int ix2_cm = (i2 - imin2) + ip2 * novp;

                        if (r1 + ri2 < rnn) {
                            // Linear interpolation
                            double xr = (rnn - ri[i1_0] - r1) / (r2 - r1);
                            int ind2 = i2 + 1;
                            if (i2 == imt[ip2]) ind2 = i2 - 1;
                            int ind2_0 = ind2 - 1;
                            xr = xr * (r2 - ri[i2_0]) / (ri[ind2_0] - ri[i2_0]);

                            cmovp[ix1_cm + istatx_dim * ix2_cm] +=
                                CFloat(static_cast<float>(temp * (1.0 - xr)), 0.0f);

                            int ix2b = (ind2 - imin2) + ip2 * novp;
                            cmovp[ix1_cm + istatx_dim * ix2b] +=
                                CFloat(static_cast<float>(temp * xr), 0.0f);

                            // r1 = rnn - ri2;  // update r1 (note: variable not used after)
                        } else {
                            cmovp[ix1_cm + istatx_dim * ix2_cm] +=
                                CFloat(static_cast<float>(temp), 0.0f);
                        }
                    } // end i2 loop (label 170)
                } // end i1 loop (label 180)
            }
        } // end iat loop (label 190)

        xn += xnatph[ip1];
    } // end ip1 loop (label 200)

    // Fill last row of cmovp using bmat
    if (ipot == 0) {
        for (int iph = 0; iph <= nph; iph++) {
            double aa = xnatph[iph] / xn;
            for (int ix1_idx = 0; ix1_idx < ncp - 1; ix1_idx++) {
                // cmovp[ncp-1 + istatx_dim * ix1_idx] is the last row
                cmovp[(ncp - 1) + istatx_dim * ix1_idx] +=
                    CFloat(static_cast<float>(aa) * bmat[iph + (nphx + 1) * ix1_idx], 0.0f);
            }
        }
    } else {
        int iph = 0;
        for (int ix1_idx = 0; ix1_idx < ncp - 1; ix1_idx++) {
            cmovp[(ncp - 1) + istatx_dim * ix1_idx] +=
                CFloat(bmat[iph + (nphx + 1) * ix1_idx], 0.0f);
        }
    }

    // LU decomposition
    int istatx = istatx_dim;
    int info = 0;
    cgetrf_(&ncp, &ncp, cmovp, &istatx, ipiv, &info);
    if (info != 0) {
        log.wlog("    *** Error in cgetrf when computing cmovp");
    }

    // Check that last pivot was not permuted
    // (ipiv is 1-based from LAPACK)
    if (ipiv[ncp - 1] != ncp) {
        feff::par::par_stop("illegal permutation in ipiv ");
    }
}

} // namespace feff::pot
