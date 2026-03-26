// Non-diagonal Lagrange parameters for SCF calculation.
// Converted from: src/ATOM/lagdat.f
//
// INDEX CONVENTION: Loop indices remain 1-based as in Fortran.
// Array accesses use F() macro for 0-based C++ arrays.

#include "lagdat.hpp"
#include "coulomb_integrals.hpp"
#include <cmath>

namespace feff::atom {

#define F(arr, i) (arr)[(i) - 1]

void lagdat(int ia, int iex, AtomState& state)
{
    auto& scf     = state.scf;
    auto& config  = state.config;
    auto& lagrange = state.lagrange;
    auto& angular = state.angular;

    int norbsc = scf.norbsc;

    // Local arrays for orbital pair indices
    int ni[2], nj[2];  // 0-based storage for 2-element arrays

    int i1 = (ia > 1) ? ia : 1;  // Fortran: i1 = max(ia, 1)
    int idep = 1;
    if (ia > 0) goto label_15;

label_11:
    idep = i1 + 1;

label_15:
    ni[0] = i1;       // ni(1) = i1
    nj[1] = i1;       // nj(2) = i1
    int ji1 = 2 * std::abs(F(config.kap, i1)) - 1;

    for (int i2 = idep; i2 <= norbsc; ++i2) {
        if (i2 == i1 || F(config.kap, i2) != F(config.kap, i1)) continue;
        if (F(lagrange.nre, i1) < 0 && F(lagrange.nre, i2) < 0) continue;
        // Skip case where occupation numbers are equal
        if (F(config.xnel, i1) == F(config.xnel, i2)) continue;

        ni[1] = i2;   // ni(2) = i2
        nj[0] = i2;   // nj(1) = i2
        double d = 0.0;

        for (int l = 1; l <= norbsc; ++l) {
            int k = 0;
            int jjl = 2 * std::abs(F(config.kap, l)) - 1;
            int kma = (ji1 < jjl) ? ji1 : jjl;  // min(ji1, jjl)

        label_41:
            {
                // akeato/bkeato/fdrirk expect 0-based orbital indices;
                // l, i1, i2 are 1-based (F() macro convention)
                double a = akeato(l - 1, i1 - 1, k, angular) / F(config.xnel, i1);
                double b = a - akeato(l - 1, i2 - 1, k, angular) / F(config.xnel, i2);
                double c_val = b;
                if (a != 0.0) c_val = c_val / a;
                if (std::abs(c_val) >= 1.0e-07) {
                    d += b * fdrirk(l - 1, l - 1, i1 - 1, i2 - 1, k, state);
                }
            }
            k = k + 2;
            if (k <= kma) goto label_41;

            if (iex == 0) continue;  // skip exchange terms

            kma = (ji1 + jjl) / 2;
            k = std::abs(jjl - kma);
            if (F(config.kap, i1) * F(config.kap, l) < 0) k = k + 1;

        label_61:
            {
                double a = bkeato(l - 1, i2 - 1, k, angular) / F(config.xnel, i2);
                double b = a - bkeato(l - 1, i1 - 1, k, angular) / F(config.xnel, i1);
                double c_val = b;
                if (a != 0.0) c_val = c_val / a;
                if (std::abs(c_val) >= 1.0e-07) {
                    d += b * fdrirk(i1 - 1, l - 1, i2 - 1, l - 1, k, state);
                }
            }
            k = k + 2;
            if (k <= kma) goto label_61;
        }

        // Store in triangular array: eps(min(i1,i2) + (max(i1,i2)-1)*(max(i1,i2)-2)/2)
        int i_min = (i1 < i2) ? i1 : i2;
        int j_max = (i1 > i2) ? i1 : i2;
        int idx = i_min + ((j_max - 1) * (j_max - 2)) / 2;
        F(lagrange.eps, idx) = d / (F(config.xnel, i2) - F(config.xnel, i1));
    }

    if (ia > 0) return;  // goto 999
    i1 = i1 + 1;
    if (i1 < norbsc) goto label_11;
}

#undef F

} // namespace feff::atom
