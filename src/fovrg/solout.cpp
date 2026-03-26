// Outward solution of the Dirac equation with power-series start.
// Converted from: src/FOVRG/solout.f
//
// INDEXING CONVENTION:
//   All array indices are 0-based. Parameters jri, max0, iwkb retain
//   their Fortran 1-based numeric values. Array element Fortran arr(k)
//   is accessed as arr[k-1] in C++.

#include "solout.hpp"
#include "intout.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <cmath>
#include <algorithm>

// Forward declaration of flatv (defined in dfovrg.cpp)
namespace feff::fovrg {
void flatv(double r1, double r2, FeffComplex p1, FeffComplex q1,
           FeffComplex en, FeffComplex vav, int ikap,
           FeffComplex& p2, FeffComplex& q2);
}

namespace feff::fovrg {

void solout(FeffComplex en, FeffComplex fl, FeffComplex agi, FeffComplex api,
            int kap, int jri, int max0, int ic3, const FeffComplex vm[],
            int iwkb, DiracWorkspaceComplex& work, MeshParamsComplex& mesh)
{
    constexpr double ccl_val = 2.0 * alpinv;
    constexpr double csq = ccl_val * ccl_val;

    double cl = work.cl;
    double hx = mesh.hx;
    int ndor = mesh.ndor;
    double* dr = mesh.dr;

    FeffComplex* gg = work.gg;
    FeffComplex* ag = work.ag;
    FeffComplex* gp = work.gp;
    FeffComplex* ap = work.ap;
    FeffComplex* dv = work.dv;
    FeffComplex* av = work.av;
    FeffComplex* ceg = work.ceg;
    FeffComplex* cep = work.cep;

    // Determine api from av[0] (Fortran av(1)) and kap
    if (av[0].real() < 0.0 && kap > 0)
        api = -agi * (FeffComplex(kap, 0.0) + fl) / av[0];
    if (av[0].real() < 0.0 && kap < 0)
        api = -agi * av[0] / (FeffComplex(kap, 0.0) - fl);

    FeffComplex ec = en / cl;
    // Fortran: ag(1) = agi, ap(1) = api => C++: ag[0], ap[0]
    ag[0] = agi;
    ap[0] = api;
    // Fortran: do i=2,ndor => ag(i) = ceg(i-1) => C++: ag[i-1] = ceg[i-2] for i=2..ndor
    for (int i = 1; i < ndor; i++) {
        ag[i] = ceg[i - 1];
        ap[i] = cep[i - 1];
    }

    if (ic3 == 0) {
        // Desclaux power expansion
        // Fortran: do j=2,ndor => C++: j=1..ndor-1 (0-based j maps to Fortran j+1)
        for (int j = 1; j < ndor; j++) {
            // Fortran k = j-1 where j is Fortran j (= C++ j + 1)
            // So Fortran k = (C++ j + 1) - 1 = C++ j
            // Fortran a = fl + kap + k = fl + kap + j
            // Fortran b = fl - kap + k = fl - kap + j
            FeffComplex a_val = fl + FeffComplex(kap + j, 0.0);
            FeffComplex b_val = fl - FeffComplex(kap, 0.0) + FeffComplex(j, 0.0);
            FeffComplex eph = a_val * b_val + av[0] * av[0];

            // Fortran: f = (ec+ccl)*ap(k) + ap(j)
            // Fortran k = C++ j => ap(k) = ap[j-1] (0-based)
            // Fortran j (= C++ j + 1) => ap(j) = ap[j] (0-based)
            FeffComplex f_val = (ec + FeffComplex(ccl_val, 0.0)) * ap[j - 1] + ap[j];
            FeffComplex g_val = ec * ag[j - 1] + ag[j];

            // Fortran: do i=1,k => f = f - av(i+1)*ap(j-i)
            // k = C++ j; Fortran i=1..k
            // av(i+1) => av[i] (C++ 0-based, i goes 1..j)
            // Fortran ap(j-i) where Fortran j = C++ j + 1, Fortran i = loop var
            // ap(Fortran_j - Fortran_i) = ap[Fortran_j - Fortran_i - 1] = ap[C++ j + 1 - (ii+1) - 1] = ap[j - 1 - ii]
            for (int ii = 0; ii < j; ii++) {
                f_val = f_val - av[ii + 1] * ap[j - 1 - ii];
                g_val = g_val - av[ii + 1] * ag[j - 1 - ii];
            }
            ag[j] = (b_val * f_val + av[0] * g_val) / eph;
            ap[j] = (av[0] * f_val - a_val * g_val) / eph;
        }

        // Initial values: Fortran do i=1,1 => gg(1) = sum
        // C++: gg[0] = sum
        gg[0] = FeffComplex(0.0, 0.0);
        gp[0] = FeffComplex(0.0, 0.0);
        for (int j = 0; j < ndor; j++) {
            // Fortran: a = fl + j - 1 (j=1..ndor) => a = fl + 0, fl+1, ...
            // C++: a = fl + j (j=0..ndor-1) => same sequence
            FeffComplex a_val = fl + FeffComplex(j, 0.0);
            // Fortran: b = dr(1)**a => C++: b = dr[0]**a
            FeffComplex b_val = std::pow(dr[0], a_val.real());
            gg[0] = gg[0] + b_val * ag[j];
            gp[0] = gp[0] + b_val * ap[j];
        }
    } else {
        // c3-corrected expansion (see fovrg.f in feff6)
        double twoz = -av[0].real() * 2.0 * cl;
        double rat1 = twoz / ccl_val;
        double rat2 = rat1 * rat1;
        double rat3 = csq / twoz;
        int il = -kap;
        if (kap > 0) il = kap + 1;
        int l0 = il - 1;

        ag[0] = agi;
        if (twoz <= 0.0) {
            ap[0] = -ec / FeffComplex(2.0 * il + 1.0, 0.0) * dr[0] * ag[0];
            ag[1] = FeffComplex(0.0, 0.0);
            ap[1] = FeffComplex(0.0, 0.0);
            ag[2] = FeffComplex(0.0, 0.0);
            ap[2] = FeffComplex(0.0, 0.0);
        } else {
            ap[0] = (fl - FeffComplex(il, 0.0)) * rat3 * ag[0];
            ag[1] = (3.0 * fl - FeffComplex(rat2, 0.0)) / (2.0 * fl + FeffComplex(1.0, 0.0)) * ag[0];
            ap[1] = rat3 * ((fl - FeffComplex(l0, 0.0)) * ag[1] - ag[0]) - ap[0];
            ag[2] = ((fl + FeffComplex(3.0 * il, 0.0)) * ag[1] - FeffComplex(3.0 * l0, 0.0) * ag[0] +
                     (fl + FeffComplex(il + 3.0, 0.0)) / rat3 * ap[1]) / (fl + FeffComplex(1.0, 0.0)) / 4.0;
            ap[2] = (rat3 * (FeffComplex(2.0 * l0, 0.0) * (fl + FeffComplex(2.0 - il, 0.0)) -
                     FeffComplex(l0 + rat2, 0.0)) * ag[1]
                     - FeffComplex(3.0 * l0, 0.0) * rat3 * (fl + FeffComplex(2.0 - il, 0.0)) * ag[0]
                     + (fl + FeffComplex(3.0 - 2.0 * il - rat2, 0.0)) * ap[1]) /
                    (fl + FeffComplex(1.0, 0.0)) / 4.0;

            ap[0] = ap[0] / ccl_val;
            ag[1] = ag[1] * rat3;
            ap[1] = ap[1] * rat3 / ccl_val;
            ag[2] = ag[2] * rat3 * rat3;
            ap[2] = ap[2] * rat3 * rat3 / ccl_val;
        }
        gg[0] = std::pow(dr[0], fl.real()) * (ag[0] + dr[0] * (ag[1] + dr[0] * ag[2]));
        gp[0] = std::pow(dr[0], fl.real()) * (ap[0] + dr[0] * (ap[1] + dr[0] * ap[2]));
    }

    // Outward integration
    // Fortran: i0 = 1 (1-based), iflat = min(jri, iwkb) (1-based)
    // intout expects i0 and max0 as 1-based values
    int i0 = 1;
    int iflat = std::min(jri, iwkb);
    intout(en, i0, kap, iflat, ic3, vm, work, mesh);

    // Use flatv beyond iflat
    // Fortran: do i = iflat, max0-1  (1-based)
    // C++ 0-based: i goes from iflat-1 to max0-2
    for (int i = iflat - 1; i <= max0 - 2; i++) {
        FeffComplex eph;
        // Fortran: if (i.eq.iwkb) => C++: if (i+1 == iwkb) since i is 0-based, iwkb is 1-based
        if (i + 1 == iwkb) {
            // Fortran: eph = cl*(3*dv(iwkb+1) - dv(iwkb+2))/2
            // dv(iwkb+1) => dv[iwkb]  (since iwkb is 1-based)
            // dv(iwkb+2) => dv[iwkb+1]
            eph = cl * (3.0 * dv[iwkb] - dv[iwkb + 1]) / 2.0;
            // Fortran: if (iwkb.eq.jri-1) => same comparison (both 1-based)
            if (iwkb == jri - 1)
                eph = cl * (dv[i] + dv[i + 1]) / 2.0;
        } else {
            eph = cl * (dv[i] + dv[i + 1]) / 2.0;
        }
        // Fortran: if (ic3.gt.0 .and. i.lt.jri)
        // i is 0-based, jri is 1-based => Fortran i < jri means C++ i+1 < jri => i < jri-1
        if (ic3 > 0 && i < jri - 1) {
            double rav = (dr[i] + dr[i + 1]) / 2.0;
            FeffComplex ec_loc = std::pow(rav, 3.0) * std::pow(FeffComplex(ccl_val, 0.0) + (en - eph) / cl, 2.0);
            eph = eph + FeffComplex(ic3, 0.0) * cl / ec_loc * (vm[i] + vm[i + 1]) / 2.0;
        }
        FeffComplex p2, q2;
        flatv(dr[i], dr[i + 1], gg[i], gp[i], en, eph, kap, p2, q2);
        gg[i + 1] = p2;
        gp[i + 1] = q2;
    }
}

} // namespace feff::fovrg
