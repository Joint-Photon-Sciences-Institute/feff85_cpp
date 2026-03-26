// Adjust hydrogen atom positions for better MT geometry.
// Converted from src/POT/moveh.f

#include "moveh.hpp"
#include "../math/distance.hpp"
#include <cmath>

namespace feff::pot {

void moveh(int nat, const int* iphat, const int* iz, double* rath)
{
    // rath is stored as rath[3][natx] in column-major (Fortran) order
    // Access: rath[j * 3 + i] for coordinate i of atom j (both 0-based)
    auto rat = [&](int coord, int iat) -> double& {
        return rath[iat * 3 + coord];
    };

    for (int iat = 0; iat < nat; ++iat) {
        if (iz[iphat[iat]] != 1) continue;

        // Find nearest atom A (units are bohr)
        double rah = 100.0;
        int ia = -1;
        for (int i = 0; i < nat; ++i) {
            if (i == iat) continue;
            double rattmp = feff::math::dist(&rath[iat * 3], &rath[i * 3]);
            if (rattmp < rah) {
                ia = i;
                rah = rattmp;
            }
        }
        if (ia < 0) continue;
        if (iz[iphat[ia]] == 1) continue;

        // Set max distance as function of rah (from H2O and GeH4 calculations)
        double ratmax = rah + 4.0 / (rah * rah);

        // Find shortest AB bond (neither A nor B is H)
        double rab = 10.0;
        int ib = -1;
        for (int i = 0; i < nat; ++i) {
            if (i == ia) continue;
            if (iz[iphat[i]] == 1) continue;
            double rattmp = feff::math::dist(&rath[ia * 3], &rath[i * 3]);
            if (rab > rattmp) {
                rab = rattmp;
                ib = i;
            }
        }
        if (rab < ratmax) ratmax = 0.95 * rab + 0.05 * rah;
        if (rah > ratmax) continue;

        // Increase rah to ratmax and check that A is still closest to H
        double ratmin = rah;
        bool retry = true;
        while (retry) {
            retry = false;
            for (int i = 0; i < 3; ++i) {
                rat(i, iat) = rat(i, ia) + ratmax / ratmin * (rat(i, iat) - rat(i, ia));
            }

            double rbh = 10.0;
            ib = -1;
            for (int i = 0; i < nat; ++i) {
                if (i == iat) continue;
                double rattmp = feff::math::dist(&rath[iat * 3], &rath[i * 3]);
                if (rbh > rattmp) {
                    rbh = rattmp;
                    ib = i;
                }
            }

            if (ia != ib) {
                rab = feff::math::dist(&rath[ia * 3], &rath[ib * 3]);
                double rattmp = ratmax * rab * rab / (ratmax * ratmax + rab * rab - rbh * rbh);
                ratmin = ratmax;
                ratmax = 0.95 * rattmp + 0.05 * rah;
                retry = true;
            }
        }
    }
}

} // namespace feff::pot
