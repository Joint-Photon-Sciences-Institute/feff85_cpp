// Angular exchange coefficients for FOVRG.
// Converted from: src/FOVRG/muatcc.f

#include "muatcc.hpp"
#include "../../src/math/wigner.hpp"
#include <feff/dimensions.hpp>
#include <cmath>

namespace feff::fovrg {

void muatcc(const double xnval[30], FovrgState& state)
{
    auto& config = state.config;
    auto& scf = state.scf;
    auto& angular = state.angular;
    int norb = scf.norb;

    // Zero the angular coefficients array
    // Fortran: do i=-ltot-1,ltot; do j=1,30; do k=0,3
    // Using the AngularCoefficientsC operator() which handles the offset
    for (int i = -ltot - 1; i <= ltot; i++) {
        for (int j = 0; j < 30; j++) {
            for (int k = 0; k < 4; k++) {
                angular(i, j, k) = 0.0;
            }
        }
    }

    // Calculate b_k(ikap, j) using Wigner 3j symbols
    for (int ikap = -ltot - 1; ikap <= ltot; ikap++) {
        if (ikap != 0) {
            int li = std::abs(ikap) * 2 - 1;
            for (int j = 0; j < norb - 1; j++) {
                int lj = std::abs(config.kap[j]) * 2 - 1;
                int kmax = (li + lj) / 2;
                int kmin = std::abs(li - lj) / 2;
                if ((ikap * config.kap[j]) < 0) kmin = kmin + 1;

                if (xnval[j] <= 0.0) {
                    // Calculate b_k(i,j)
                    for (int k = kmin; k <= kmax; k += 2) {
                        int index = (k - kmin) / 2;
                        double w3j = feff::math::cwig3j(li, k * 2, lj, 1, 0, 2);
                        angular(ikap, j, index) = config.xnel[j] * (w3j * w3j);
                    }
                }
            }
        }
    }
}

} // namespace feff::fovrg
