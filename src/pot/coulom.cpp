// Coulomb potential with charge normalization.
// Converted from src/POT/coulom.f

#include "coulom.hpp"
#include "frnrm.hpp"
#include "../atom/utility.hpp"
#include "../math/distance.hpp"
#include <feff/constants.hpp>
#include <cmath>

namespace feff::pot {

double fab(double aa, double bb, double r0, double r1, double r2)
{
    double a2 = (r2 * r2 - r1 * r1) / 2.0;
    double a3 = (r2 * r2 * r2 - r1 * r1 * r1) / 3.0;
    double a4 = (r2 * r2 * r2 * r2 - r1 * r1 * r1 * r1) / 4.0;
    return aa * (a4 / r0 - a3) + bb * (a3 / r0 - a2);
}

void coulom(int icoul, int npot, const int* ilast, const double* rhoval,
            const double* edenvl, const double* edens,
            int nat, const double* rat, const int* iatph, const int* iphat,
            double* rnrm, const double* dq, const int* iz, double* vclap)
{
    constexpr int s251 = 251;
    double ri05[s251];
    double drho[s251], dvcl[s251];

    // Make radial grid with 0.05 step
    double dx05 = 0.05;
    for (int i = 0; i < s251; ++i) {
        ri05[i] = std::exp(-8.8 + dx05 * i);
    }

    for (int ip = 0; ip <= npot; ++ip) {
        for (int ir = 0; ir < ilast[ip]; ++ir) {
            drho[ir] = (rhoval[ip * s251 + ir] - edenvl[ip * s251 + ir]) *
                       ri05[ir] * ri05[ir];
        }

        feff::atom::potslw(dvcl, drho, ri05, dx05, ilast[ip]);

        for (int ir = ilast[ip]; ir < s251; ++ir) {
            dvcl[ir] = 0.0;
        }

        double dvnrm;
        if (icoul == 1) {
            // Explicit charge normalization
            int jnrm = static_cast<int>((std::log(rnrm[ip]) + 8.8) / 0.05) + 2;
            dvnrm = dq[ip] / rnrm[ip];

            int iat0 = iatph[ip];  // 1-based
            for (int iat = 0; iat < nat; ++iat) {
                if (iat == iat0 - 1) continue;
                double rr_dist = feff::math::dist(&rat[iat * 3], &rat[(iat0 - 1) * 3]);
                if (rr_dist < rnrm[ip]) rr_dist = rnrm[ip];
                dvnrm += dq[iphat[iat]] / rr_dist;
            }

            // Transfer condition to r(jnrm) instead of rnrm
            double dr = ri05[jnrm - 1] - rnrm[ip];
            double bb = (drho[jnrm - 1] - drho[jnrm - 2]) /
                        (ri05[jnrm - 1] - ri05[jnrm - 2]);
            dvnrm = dvnrm - dr / 2.0 * (dq[ip] / (rnrm[ip] * rnrm[ip]) +
                    (dq[ip] + drho[jnrm - 1] * dr - bb / 2.0 * dr * dr) /
                    (ri05[jnrm - 1] * ri05[jnrm - 1]));

            dvnrm = dvnrm - dvcl[jnrm - 1];
        } else {
            // Norman picture normalization (default, icoul=0)
            double rnrm1;
            frnrm(&edens[ip * s251], iz[ip], rnrm1);

            double drho_full[s251];
            for (int i = 0; i < s251; ++i) {
                drho_full[i] = edens[ip * s251 + i] - edenvl[ip * s251 + i] +
                               rhoval[ip * s251 + i];
            }
            double rnrm2;
            frnrm(drho_full, iz[ip], rnrm2);

            double rmin = std::min(rnrm1, rnrm2);
            int inrm = static_cast<int>((std::log(rmin) + 8.8) / 0.05) + 1;
            double r0 = ri05[inrm - 1];

            double delv = 0.0;
            if (rnrm2 > rnrm1) {
                double aa = (drho_full[inrm] - drho_full[inrm - 1]) /
                            (ri05[inrm] - ri05[inrm - 1]);
                double bb = drho_full[inrm - 1] - aa * ri05[inrm - 1];
                delv -= fab(aa, bb, r0, rnrm1, rnrm2);
            } else {
                double aa = (edens[ip * s251 + inrm - 1] - edens[ip * s251 + inrm]) /
                            (ri05[inrm] - ri05[inrm - 1]);
                double bb = -edens[ip * s251 + inrm - 1] - aa * ri05[inrm - 1];
                delv -= fab(aa, bb, r0, rnrm2, rnrm1);
            }

            double aa = (drho_full[inrm] - drho_full[inrm - 1] +
                         edens[ip * s251 + inrm - 1] - edens[ip * s251 + inrm]) /
                        (ri05[inrm] - ri05[inrm - 1]);
            double bb = drho_full[inrm - 1] - edens[ip * s251 + inrm - 1] -
                        aa * ri05[inrm - 1];
            delv -= fab(aa, bb, r0, r0, rmin);

            dvnrm = delv - dvcl[inrm - 1];
        }

        for (int ir = 0; ir < ilast[ip]; ++ir) {
            vclap[ip * s251 + ir] += dvcl[ir] + dvnrm;
        }
        for (int ir = ilast[ip]; ir < s251; ++ir) {
            vclap[ip * s251 + ir] = 0.0;
        }
    }
}

} // namespace feff::pot
