// Outward integration of the Dirac equation (complex version).
// Converted from: src/FOVRG/intout.f
//
// Runge-Kutta start for first npi points, then Milne predictor-corrector.
//
// INDEXING CONVENTION:
//   All array indices are 0-based (C++ standard).
//   Parameters i0 and max0 are 1-based (Fortran convention) and are
//   converted to 0-based at the top of this function.
//   The grid arrays dr[], dv[], gg[], gp[], eg[], ep[], vm[] use 0-based
//   indexing where C++ index k corresponds to Fortran index k+1.

#include "intout.hpp"
#include <feff/dimensions.hpp>
#include <cmath>

namespace feff::fovrg {

void intout(FeffComplex en, int i0, int kap, int max0, int ic3,
            const FeffComplex vm[], DiracWorkspaceComplex& work,
            MeshParamsComplex& mesh)
{
    constexpr int npi = 6;
    constexpr double test = 1.0e+5;

    double cl = work.cl;
    double hx = mesh.hx;
    int np = mesh.np;
    double* dr = mesh.dr;

    // Aliases for workspace arrays
    FeffComplex* gg = work.gg;
    FeffComplex* ag = work.ag;
    FeffComplex* gp = work.gp;
    FeffComplex* ap = work.ap;
    FeffComplex* dv = work.dv;
    FeffComplex* av = work.av;
    FeffComplex* eg = work.eg;
    FeffComplex* ep = work.ep;

    // i0 and max0 come in as 1-based Fortran values.
    // Convert to 0-based for all array access.
    // Fortran: i = i0, uses dv(i), dr(i)  => C++: i = i0-1, uses dv[i], dr[i]
    // Fortran: max0 is last point to integrate to
    // Fortran loop "do i = ..., max0-1" computes gg(max0); C++ computes gg[max0-1].

    double ccl = cl + cl;
    double exphx = std::exp(hx / 2.0);
    int ihard = 0;
    FeffComplex ec = en / cl;

    // Local derivative arrays for Runge-Kutta / Milne (1-based in Fortran, 0-based here)
    FeffComplex dg_loc[npi], dp_loc[npi];

    // Runge-Kutta for first npi points
    // Fortran: i = i0, j = 1
    // C++: i = i0 - 1 (0-based), j = 0 (0-based into dg_loc)
    int i = i0 - 1;
    int j = 0;

    FeffComplex f = (ec - dv[i]) * dr[i];
    FeffComplex g = f + ccl * dr[i];
    FeffComplex c3 = FeffComplex(ic3, 0.0) * vm[i] / (g * g);
    dg_loc[j] = hx * (g * gp[i] - FeffComplex(kap, 0.0) * gg[i] + ep[i]);
    dp_loc[j] = hx * (FeffComplex(kap, 0.0) * gp[i] - (f - c3) * gg[i] - eg[i]);

    // Runge-Kutta loop
    // Fortran: if (i.lt.max0) then ...
    // max0 is 1-based; i is 0-based. Fortran i < max0 means i_0based < max0 - 1.
    while (i < max0 - 1) {
        FeffComplex ac = gg[i] + 0.5 * dg_loc[j];
        FeffComplex bc = gp[i] + 0.5 * dp_loc[j];
        double rh = dr[i] * exphx;

        // Interpolate potential and exchange between two points
        // Fortran: dr(i+1) => C++: dr[i+1] (next grid point)
        double xm1 = (dr[i + 1] - rh) / (dr[i + 1] - dr[i]);
        double xm2 = (rh - dr[i]) / (dr[i + 1] - dr[i]);

        FeffComplex vh, vmh;
        // Fortran: i0.eq.1 => C++: i0 == 1 (i0 is still the original 1-based value)
        if (av[0].real() < 0.0 && i0 == 1) {
            // Point nucleus: important nonlinearity from z/r term
            FeffComplex dv1 = dv[i] - av[0] / dr[i];
            FeffComplex dv2 = dv[i + 1] - av[0] / dr[i + 1];
            vh = dv1 * xm1 + dv2 * xm2;
            vh = vh + av[0] / rh;
            vmh = (xm1 * vm[i] * dr[i] + xm2 * vm[i + 1] * dr[i + 1]) / rh;
        } else if (i0 == 1) {
            // Finite nucleus: nonlinearity from z*r^2 term
            // Fortran: av(4) => C++: av[3]
            FeffComplex dv1 = dv[i] - av[3] * dr[i] * dr[i];
            FeffComplex dv2 = dv[i + 1] - av[3] * dr[i + 1] * dr[i + 1];
            vh = (dv1 * (dr[i + 1] - rh) + dv2 * (rh - dr[i])) / (dr[i + 1] - dr[i]);
            vh = vh + av[3] * rh * rh;
            vmh = (xm1 * vm[i] / (dr[i] * dr[i]) + xm2 * vm[i + 1] / (dr[i + 1] * dr[i + 1])) * rh * rh;
        } else {
            // Outward integration of irregular solution near jri
            vh = dv[i] * xm1 + dv[i + 1] * xm2;
            vmh = xm1 * vm[i] + xm2 * vm[i + 1];
        }
        FeffComplex eph = ep[i] * xm1 + ep[i + 1] * xm2;
        FeffComplex egh = eg[i] * xm1 + eg[i + 1] * xm2;

        f = (ec - vh) * rh;
        g = f + ccl * rh;
        c3 = FeffComplex(ic3, 0.0) * vmh / (g * g);
        FeffComplex dg2 = hx * (g * bc - FeffComplex(kap, 0.0) * ac + eph);
        FeffComplex dp2 = hx * (FeffComplex(kap, 0.0) * bc - (f - c3) * ac - egh);
        ac = ac + 0.5 * (dg2 - dg_loc[j]);
        bc = bc + 0.5 * (dp2 - dp_loc[j]);
        FeffComplex dg3 = hx * (g * bc - FeffComplex(kap, 0.0) * ac + eph);
        FeffComplex dp3 = hx * (FeffComplex(kap, 0.0) * bc - (f - c3) * ac - egh);
        ac = ac + dg3 - 0.5 * dg2;
        bc = bc + dp3 - 0.5 * dp2;

        i = i + 1;
        j = j + 1;
        f = (ec - dv[i]) * dr[i];
        g = f + ccl * dr[i];
        c3 = FeffComplex(ic3, 0.0) * vm[i] / (g * g);
        FeffComplex dg4 = hx * (g * bc - FeffComplex(kap, 0.0) * ac + ep[i]);
        FeffComplex dp4 = hx * (FeffComplex(kap, 0.0) * bc - (f - c3) * ac - eg[i]);
        gg[i] = gg[i - 1] + (dg_loc[j - 1] + 2.0 * (dg2 + dg3) + dg4) / 6.0;
        gp[i] = gp[i - 1] + (dp_loc[j - 1] + 2.0 * (dp2 + dp3) + dp4) / 6.0;
        dg_loc[j] = hx * (g * gp[i] - FeffComplex(kap, 0.0) * gg[i] + ep[i]);
        dp_loc[j] = hx * (FeffComplex(kap, 0.0) * gp[i] - (f - c3) * gg[i] - eg[i]);
        if (j < npi - 1) continue;

        // Scale derivatives for Milne method
        for (int ii = 0; ii < npi; ii++) {
            dg_loc[ii] = dg_loc[ii] / hx;
            dp_loc[ii] = dp_loc[ii] / hx;
        }

        // Integration of the inhomogeneous system by Milne predictor-corrector
        double a1 = hx * 3.3;
        double a2 = -hx * 4.2;
        double a3 = hx * 7.8;
        double a4 = hx * 14.0 / 45.0;
        double a5 = hx * 64.0 / 45.0;
        double a6 = hx * 24.0 / 45.0;

        // Fortran: do i = npi+i0-1, max0-1
        // i0 is 1-based; npi+i0-1 = npi+1-1 = npi (1-based Fortran index)
        // 0-based: npi - 1
        // max0-1 is 1-based => 0-based: max0 - 2
        // Fortran computes gg(i+1) in body, so last computed = gg(max0), 0-based: gg[max0-1]
        for (int ii = i0 - 1 + npi - 1; ii <= max0 - 2; ii++) {
            int nit = 0;
            // Predictor: Fortran gg(i-5) => C++ gg[ii-5]
            FeffComplex acp = gg[ii - 5] + a1 * (dg_loc[npi - 1] + dg_loc[npi - 5]) +
                              a2 * (dg_loc[npi - 2] + dg_loc[npi - 4]) + a3 * dg_loc[npi - 3];
            FeffComplex bcp = gp[ii - 5] + a1 * (dp_loc[npi - 1] + dp_loc[npi - 5]) +
                              a2 * (dp_loc[npi - 2] + dp_loc[npi - 4]) + a3 * dp_loc[npi - 3];

            // Corrector (without contribution from i+1)
            FeffComplex ac_m = gg[ii - 3] + a4 * dg_loc[npi - 4] +
                               a5 * (dg_loc[npi - 1] + dg_loc[npi - 3]) + a6 * dg_loc[npi - 2];
            FeffComplex bc_m = gp[ii - 3] + a4 * dp_loc[npi - 4] +
                               a5 * (dp_loc[npi - 1] + dp_loc[npi - 3]) + a6 * dp_loc[npi - 2];

            // Shift derivative arrays
            for (int jj = 0; jj < npi - 1; jj++) {
                dg_loc[jj] = dg_loc[jj + 1];
                dp_loc[jj] = dp_loc[jj + 1];
            }

            f = (ec - dv[ii + 1]) * dr[ii + 1];
            g = f + ccl * dr[ii + 1];
            c3 = FeffComplex(ic3, 0.0) * vm[ii + 1] / (g * g);

            // Iterate corrector
            while (true) {
                dg_loc[npi - 1] = g * bcp - FeffComplex(kap, 0.0) * acp + ep[ii + 1];
                dp_loc[npi - 1] = FeffComplex(kap, 0.0) * bcp - (f - c3) * acp - eg[ii + 1];

                gg[ii + 1] = ac_m + a4 * dg_loc[npi - 1];
                gp[ii + 1] = bc_m + a4 * dp_loc[npi - 1];

                if (std::abs(test * (gg[ii + 1] - acp)) > std::abs(gg[ii + 1]) ||
                    std::abs(test * (gp[ii + 1] - bcp)) > std::abs(gp[ii + 1])) {
                    if (nit < 40) {
                        acp = gg[ii + 1];
                        bcp = gp[ii + 1];
                        nit++;
                        continue;
                    } else {
                        ihard++;
                    }
                }
                break;
            }
        }
        break;  // Exit the RK loop after Milne section completes
    }

    // Zero out beyond max0
    // Fortran: do i = max0+1, np => zeroes gg(max0+1)..gg(np)
    // 0-based: gg[max0]..gg[np-1]
    for (int ii = max0; ii < np; ii++) {
        gg[ii] = FeffComplex(0.0, 0.0);
        gp[ii] = FeffComplex(0.0, 0.0);
    }
}

} // namespace feff::fovrg
