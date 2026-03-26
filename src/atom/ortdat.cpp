// Schmidt orthogonalization of orbitals.
// Converted from: src/ATOM/ortdat.f
//
// INDEX CONVENTION: Loop indices remain 1-based as in Fortran.
// Array accesses use F() macro for 0-based C++ arrays.
//
// Uses dsordf() for computing overlap and norm integrals.
// In the Fortran, dsordf uses the working arrays dg/ag/dp/ap from COMMON /comdir/
// to exchange data with ortdat (ortdat loads them, dsordf reads them).

#include "ortdat.hpp"
#include "dsordf.hpp"
#include <cmath>
#include <algorithm>

namespace feff::atom {

#define F(arr, i) (arr)[(i) - 1]

void ortdat(int ia, AtomState& state)
{
    auto& orb    = state.orb;
    auto& work   = state.work;
    auto& config = state.config;
    auto& scf    = state.scf;
    auto& mesh   = state.mesh;

    int norb = scf.norb;
    int ndor = mesh.ndor;
    int idim = mesh.idim;

    // Aliases to workspace arrays used for exchange with dsordf
    double* dg_work = work.dg;   // large component workspace
    double* ag_work = work.ag;
    double* dp_work = work.dp;   // small component workspace
    double* ap_work = work.ap;

    int m = norb;
    int l = (ia > 1) ? ia : 1;  // Fortran: l = max(ia, 1)

    if (ia <= 0) goto label_5;
    goto label_11;

label_5:
    m = l;
    l = l + 1;
    if (l > norb) return;  // goto 999

label_11:
    // Zero out workspace arrays
    for (int i = 1; i <= idim; ++i) {
        F(dg_work, i) = 0.0;
        F(dp_work, i) = 0.0;
    }

    {
        int maxl = F(config.nmax, l);

        // Load orbital l into workspace
        for (int i = 1; i <= maxl; ++i) {
            // cg(i,l) in Fortran -> orb.cg[i-1][l-1] in C++
            F(dg_work, i) = orb.cg[i - 1][l - 1];
            F(dp_work, i) = orb.cp[i - 1][l - 1];
        }
        for (int i = 1; i <= ndor; ++i) {
            F(ag_work, i) = orb.bg[i - 1][l - 1];
            F(ap_work, i) = orb.bp[i - 1][l - 1];
        }

        // Subtract projection onto each orbital j with same kappa
        for (int j = 1; j <= m; ++j) {
            if (j != l && F(config.kap, j) == F(config.kap, l)) {
                int max0 = F(config.nmax, j);
                // dsordf(j, j, 0, 3, fl(l)) computes overlap <j|workspace>
                // jnd=3: hg(l)=dg(l)*cg(l,i)+dp(l)*cp(l,j)
                double a = dsordf(j, j, 0, 3, orb.fl[l - 1], state);
                for (int i = 1; i <= max0; ++i) {
                    F(dg_work, i) -= a * orb.cg[i - 1][j - 1];
                    F(dp_work, i) -= a * orb.cp[i - 1][j - 1];
                }
                for (int i = 1; i <= ndor; ++i) {
                    F(ag_work, i) -= a * orb.bg[i - 1][j - 1];
                    F(ap_work, i) -= a * orb.bp[i - 1][j - 1];
                    maxl = std::max(maxl, max0);
                }
            }
        }

        int max0 = maxl;
        F(config.nmax, l) = max0;

        // Normalize: dsordf(l, max0, 0, 4, fl(l)) computes <workspace|workspace>
        double a = dsordf(l, max0, 0, 4, orb.fl[l - 1], state);
        a = std::sqrt(a);

        for (int i = 1; i <= max0; ++i) {
            orb.cg[i - 1][l - 1] = F(dg_work, i) / a;
            orb.cp[i - 1][l - 1] = F(dp_work, i) / a;
        }
        for (int i = 1; i <= ndor; ++i) {
            orb.bg[i - 1][l - 1] = F(ag_work, i) / a;
            orb.bp[i - 1][l - 1] = F(ap_work, i) / a;
        }
    }

    if (ia <= 0) goto label_5;
    // label_999: return
}

#undef F

} // namespace feff::atom
