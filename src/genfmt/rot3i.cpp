// Rotation matrix calculation.
// Converted from GENFMT/rot3i.f
//
// Calculates the beta-dependence of rotation matrix elements using
// recursion of formula (4.4.1) in Edmonds.
//
// Notation: dri0[l][m][n] = D^{l'}_{m',n'} where
//   l' = l (0-based maps to Fortran l-1),  n' = n-l, m' = m-l
// Thus dri0[0][0][0] is l'=0, m'=0, n'=0.
//      dri0[2][4][4] is l'=2, m'=2, n'=2.

#include "rot3i.hpp"
#include <cmath>
#include <algorithm>

namespace feff::genfmt {

void rot3i(int lxp1, int mxp1, int ileg, const double beta[],
           RotationMatrixData& rm, int istore) {
    // ileg: index into beta[] (1-based from rdpath)
    // istore: index for rm.dri[...][...][...][istore] storage (0-based for fmtrxi)
    // If istore < 0, use ileg directly (backward compat)

    // dri0 uses full (2l+1) x (2l+1) blocks per l.
    // Indexed 0-based throughout: dri0[il][in][im] where
    // il = 0..ltot, in/im = 0..2*ltot
    // This corresponds to Fortran dri0(il+1, in+1, im+1).
    double dri0[ltot + 1][2 * ltot + 1][2 * ltot + 1];

    // Initialize dri0
    for (int il = 0; il <= ltot; ++il)
        for (int in = 0; in < 2 * ltot + 1; ++in)
            for (int im = 0; im < 2 * ltot + 1; ++im)
                dri0[il][in][im] = 0.0;

    int nm = mxp1;
    int ndm = lxp1 + nm - 1;
    double xc = std::cos(beta[ileg] / 2.0);
    double xs = std::sin(beta[ileg] / 2.0);
    double s = std::sin(beta[ileg]);

    // l'=0: dri0(1,1,1) in Fortran -> dri0[0][0][0] in C++
    dri0[0][0][0] = 1.0;

    // l'=1: dri0(2,...) in Fortran -> dri0[1][...][...] in C++
    dri0[1][0][0] = xc * xc;
    dri0[1][0][1] = s / std::sqrt(2.0);
    dri0[1][0][2] = xs * xs;
    dri0[1][1][0] = -dri0[1][0][1];
    dri0[1][1][1] = std::cos(beta[ileg]);
    dri0[1][1][2] = dri0[1][0][1];
    dri0[1][2][0] = dri0[1][0][2];
    dri0[1][2][1] = -dri0[1][1][2];
    dri0[1][2][2] = dri0[1][0][0];

    // Recursion for l'>=2
    // Fortran: do l = 3, lxp1   (l is 1-based index, l' = l-1 >= 2)
    // C++: fl = 3..lxp1 (Fortran 1-based), C++ 0-based index = fl-1
    for (int fl = 3; fl <= lxp1; ++fl) {
        int il = fl - 1;  // 0-based C++ index

        int ln = std::min(2 * fl - 1, ndm);
        int lm = std::min(2 * fl - 3, ndm);

        // Fortran loops: n = 1..ln, m = 1..lm (all 1-based)
        // C++ 0-based: fn = 1..ln, fm = 1..lm
        for (int fn = 1; fn <= ln; ++fn) {
            for (int fm = 1; fm <= lm; ++fm) {
                double t1 = static_cast<double>((2 * fl - 1 - fn) * (2 * fl - 2 - fn));
                double t  = static_cast<double>((2 * fl - 1 - fm) * (2 * fl - 2 - fm));
                double f1 = std::sqrt(t1 / t);
                double f2 = std::sqrt(static_cast<double>((2 * fl - 1 - fn) * (fn - 1)) / t);
                double t3 = static_cast<double>((fn - 2) * (fn - 1));
                double f3 = std::sqrt(t3 / t);

                // Fortran: dri0(l, n, m) uses 1-based indexing
                // C++ 0-based: dri0[fl-1][fn-1][fm-1], dri0[fl-2][...][...]
                double dlnm = f1 * xc * xc * dri0[il - 1][fn - 1][fm - 1];
                if (fn - 1 > 0)
                    dlnm -= f2 * s * dri0[il - 1][fn - 2][fm - 1];
                if (fn - 2 > 0)
                    dlnm += f3 * xs * xs * dri0[il - 1][fn - 3][fm - 1];

                dri0[il][fn - 1][fm - 1] = dlnm;

                // Symmetry: if n > (2*l-3), dri0(l,m,n) = (-1)^(n-m) * dri0(l,n,m)
                if (fn > (2 * fl - 3)) {
                    int sign = ((fn - fm) % 2 == 0) ? 1 : -1;
                    dri0[il][fm - 1][fn - 1] = sign * dri0[il][fn - 1][fm - 1];
                }
            }

            // Special corner elements when n > (2*l-3)
            if (fn > (2 * fl - 3)) {
                // Fortran 1-based indices: (l, 2l-2, 2l-2), (l, 2l-1, 2l-2), etc.
                // C++ 0-based: [il][2*fl-3][2*fl-3], etc.
                dri0[il][2*fl - 3][2*fl - 3] = dri0[il][1][1];
                dri0[il][2*fl - 2][2*fl - 3] = -dri0[il][0][1];
                dri0[il][2*fl - 3][2*fl - 2] = -dri0[il][1][0];
                dri0[il][2*fl - 2][2*fl - 2] = dri0[il][0][0];
            }
        }
    }

    // Copy result into rm.dri[...][...][...][ileg]
    // Zero the target slice first
    for (int il = 0; il <= ltot; ++il)
        for (int m1 = 0; m1 < 2 * mtot + 1; ++m1)
            for (int m2 = 0; m2 < 2 * mtot + 1; ++m2)
                rm.dri[il][m1][m2][istore] = 0.0;

    // Fortran: dri(il, m1+mtot+1, m2+mtot+1, ileg) = dri0(il, m1+il, m2+il)
    // il is 1-based in Fortran (1..lxp1); m1,m2 range -mx..mx where mx = min(il-1, mxp1-1)
    // C++ 0-based: rm.dri[il-1][m1+mtot][m2+mtot][ileg] = dri0[il-1][m1+(il-1)][m2+(il-1)]
    for (int fl = 1; fl <= lxp1; ++fl) {
        int il = fl - 1;  // 0-based
        int mx = std::min(fl - 1, mxp1 - 1);
        for (int m1 = -mx; m1 <= mx; ++m1) {
            for (int m2 = -mx; m2 <= mx; ++m2) {
                rm.dri[il][m1 + mtot][m2 + mtot][istore] =
                    dri0[il][m1 + fl - 1][m2 + fl - 1];
            }
        }
    }
}

} // namespace feff::genfmt
