// Broyden density mixing algorithm.
// Converted from src/POT/broydn.f

#include "broydn.hpp"
#include "../math/sommerfeld.hpp"
#include <feff/constants.hpp>
#include <cmath>
#include <cstring>

namespace feff::pot {

// Maximum number of Broyden iterations stored
static constexpr int nbr = 30;
static constexpr int s251 = 251;

// Static workspace for Broyden algorithm (saved between calls)
static double cmi[nbr * nbr];
static double frho[s251 * (nphx + 1) * nbr];
static double urho[s251 * (nphx + 1) * nbr];
static double xnorm_arr[nbr];
static double wt[s251];
static double rhoold[s251 * (nphx + 1)];
static double ri05_static[s251];

void broydn(int iscmt, double ca, int nph, const double* xnvmu,
            const int* ilast, const double* xnatph, const double* rnrm,
            double* qnrm, double* edenvl, double* rhoval, double* dq)
{
    double xpc[s251];
    double dx05 = 0.05;

    // Make radial grid with 0.05 step (first call only)
    if (iscmt == 1) {
        for (int i = 0; i < s251; ++i) {
            ri05_static[i] = std::exp(-8.8 + dx05 * i);
            wt[i] = ri05_static[i] * ri05_static[i] * ri05_static[i];
        }
    }

    // iscmt is 1-based; Fortran frho(ir, iph, iscmt)
    // frho index: (iscmt-1) * s251*(nphx+1) + iph*s251 + ir
    int scmt_off = (iscmt - 1) * s251 * (nphx + 1);

    // Record F(rho_i)
    for (int iph = 0; iph <= nph; ++iph) {
        for (int ir = 0; ir < ilast[iph]; ++ir) {
            frho[scmt_off + iph * s251 + ir] =
                rhoval[iph * s251 + ir] * ri05_static[ir] -
                edenvl[iph * s251 + ir] * wt[ir];
        }
    }

    // dq = total number of valence electrons for each potential
    double xnferm = 0.0;
    for (int ip = 0; ip <= nph; ++ip) {
        dq[ip] = 0.0;
        for (int il = 0; il <= lx; ++il) {
            dq[ip] += xnvmu[ip * (lx + 1) + il];
        }
        xnferm += dq[ip] * xnatph[ip];
    }

    if (iscmt > 1) {
        int prev_off = (iscmt - 2) * s251 * (nphx + 1);

        // Get normalization factor
        xnorm_arr[iscmt - 1] = 0.0;
        for (int iph = 0; iph <= nph; ++iph) {
            for (int ir = 0; ir < ilast[iph]; ++ir) {
                double diff = frho[scmt_off + iph * s251 + ir] -
                              frho[prev_off + iph * s251 + ir];
                xnorm_arr[iscmt - 1] += diff * diff;
            }
        }

        // Calculate c_m,i
        for (int j = 1; j < iscmt; ++j) {
            int j_off = j * s251 * (nphx + 1);
            int jm1_off = (j - 1) * s251 * (nphx + 1);
            // cmi[iscmt-1][j] => cmi[(iscmt-1)*nbr + j]
            cmi[(iscmt - 1) * nbr + j] = 0.0;
            for (int iph = 0; iph <= nph; ++iph) {
                for (int ir = 0; ir < ilast[iph]; ++ir) {
                    cmi[(iscmt - 1) * nbr + j] +=
                        frho[scmt_off + iph * s251 + ir] *
                        (frho[j_off + iph * s251 + ir] - frho[jm1_off + iph * s251 + ir]);
                }
            }
            cmi[(iscmt - 1) * nbr + j] /= xnorm_arr[j];
        }

        // Calculate U_i - vector of Lagrange multipliers
        int urho_off = (iscmt - 1) * s251 * (nphx + 1);
        for (int iph = 0; iph <= nph; ++iph) {
            for (int ir = 0; ir < ilast[iph]; ++ir) {
                urho[urho_off + iph * s251 + ir] =
                    ca * (frho[scmt_off + iph * s251 + ir] -
                          frho[prev_off + iph * s251 + ir]) +
                    (edenvl[iph * s251 + ir] - rhoold[iph * s251 + ir]) * wt[ir];
            }
        }

        for (int j = 1; j < iscmt - 1; ++j) {
            int uj_off = j * s251 * (nphx + 1);
            double cm_diff = cmi[(iscmt - 1) * nbr + j] -
                             cmi[(iscmt - 2) * nbr + j];
            for (int iph = 0; iph <= nph; ++iph) {
                for (int ir = 0; ir < ilast[iph]; ++ir) {
                    urho[urho_off + iph * s251 + ir] -=
                        urho[uj_off + iph * s251 + ir] * cm_diff;
                }
            }
        }
    }

    // Construct new density
    for (int iph = 0; iph <= nph; ++iph) {
        for (int ir = 0; ir < ilast[iph]; ++ir) {
            rhoold[iph * s251 + ir] = edenvl[iph * s251 + ir];
            rhoval[iph * s251 + ir] = edenvl[iph * s251 + ir] +
                ca * frho[scmt_off + iph * s251 + ir] / wt[ir];
            for (int j = 1; j < iscmt; ++j) {
                int uj_off = j * s251 * (nphx + 1);
                rhoval[iph * s251 + ir] -=
                    cmi[(iscmt - 1) * nbr + j] *
                    urho[uj_off + iph * s251 + ir] / wt[ir];
            }
        }
    }

    // Calculate charge inside Norman sphere
    double x0 = 8.8;
    double dqav = 0.0;
    double xnat = 0.0;
    for (int iph = 0; iph <= nph; ++iph) {
        int jnrm = static_cast<int>((std::log(rnrm[iph]) + x0) / dx05) + 2;
        int i0 = jnrm + 1;
        double xirf = 2.0;
        for (int ir = 0; ir < ilast[iph]; ++ir) {
            xpc[ir] = rhoval[iph * s251 + ir] * ri05_static[ir] * ri05_static[ir];
        }
        feff::math::somm2(ri05_static, xpc, dx05, xirf, rnrm[iph], 0, i0);
        dq[iph] = xirf - qnrm[iph] - dq[iph];
        dqav += xnatph[iph] * dq[iph];
        xnat += xnatph[iph];
    }

    // Keep charge neutrality
    double aa = dqav / xnferm;
    dqav = dqav / xnat;
    for (int iph = 0; iph <= nph; ++iph) {
        dq[iph] = dq[iph] - dqav;
        qnrm[iph] = qnrm[iph] + dq[iph];
        for (int ir = 0; ir < ilast[iph]; ++ir) {
            rhoval[iph * s251 + ir] -= aa * edenvl[iph * s251 + ir];
        }
    }
}

} // namespace feff::pot
