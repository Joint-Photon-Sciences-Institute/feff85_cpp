// Angular Coulomb coefficients for the ATOM module.
// Converted from: muatco.f

#include "muatco.hpp"
#include "../math/wigner.hpp"
#include <cmath>

namespace feff::atom {

void muatco(const double xnval[], ScfParams& scf,
            AngularCoefficients& ang, OrbitalConfig& config) {
    // Zero out angular coefficients
    for (int i = 0; i < max_orb; ++i) {
        for (int j = 0; j < max_orb; ++j) {
            for (int k = 0; k < 4; ++k) {
                ang.afgk[i][j][k] = 0.0;
            }
        }
    }

    // Fortran: do i=1,norb; do j=1,i  (1-based)
    // C++: i from 0 to norb-1; j from 0 to i
    for (int i = 0; i < scf.norb; ++i) {
        int li = std::abs(config.kap[i]) * 2 - 1;

        for (int j = 0; j <= i; ++j) {
            int lj = std::abs(config.kap[j]) * 2 - 1;
            int kmax = (li + lj) / 2;
            int kmin = std::abs(li - lj) / 2;

            if (config.kap[i] * config.kap[j] < 0) kmin = kmin + 1;

            // Calculate a_k(i,j)
            int m = 0;
            if (j == i && xnval[i] <= 0.0) m = 1;

            // Fortran: afgk(j,i,0) += xnel(i)*(xnel(j)-m)
            // In Fortran afgk is (30,30,0:3) with 1-based orbital indices.
            // In C++ afgk is [max_orb][max_orb][4] with 0-based.
            // afgk(j_f, i_f, 0) where j_f <= i_f corresponds to afgk[j][i][0].
            ang.afgk[j][i][0] += config.xnel[i] * (config.xnel[j] - m);

            if (xnval[i] <= 0.0 || xnval[j] <= 0.0) {
                // Calculate b_k(i,j)
                double b_val = ang.afgk[j][i][0];

                if (j == i && xnval[i] <= 0.0) {
                    double a_val = static_cast<double>(li);
                    b_val = -b_val * (a_val + 1.0) / a_val;
                    kmin = kmin + 2;
                }

                for (int k = kmin; k <= kmax; k += 2) {
                    // Fortran: afgk(i,j,k/2) = b * cwig3j(li,k*2,lj,1,0,2)^2
                    // Note: afgk(i_f, j_f, k/2) where i_f >= j_f (i.e. max,min)
                    // stores the exchange coefficient.
                    double w = feff::math::cwig3j(li, k * 2, lj, 1, 0, 2);
                    ang.afgk[i][j][k / 2] = b_val * (w * w);
                }
            }
        }
    }
}

} // namespace feff::atom
