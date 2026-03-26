// Singularity finder for Hedin-Lundqvist self-energy integrands
// Converted from src/EXCH/fndsng.f

#include "fndsng.hpp"

#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>

#include <feff/types.hpp>
#include "../math/polynomial_roots.hpp"
#include "../common/qsort.hpp"

namespace feff::exch {

void fndsng(FeffComplex limit1, FeffComplex limit2,
            int& nsing, FeffComplex xsing[20],
            const double dppar[10], const FeffComplex cpar[10], int ifcn)
{
    // Local variables
    FeffComplex coef[4], sol[4], xsing2[4];
    int nsol;
    constexpr double zero_tol = 1.0e-4;

    nsing = 0;

    // Solve eq 1 for q with + sign:
    //   +k*q^3 + 2*(3*k^2 - E - 2/3)*q^2 + 4*k*(k^2 - E)*q + [(k^2 - E)^2 - Wp^2] = 0
    coef[0] = 4.0 * cpar[0];
    coef[1] = 2.0 * (3.0 * cpar[0] * cpar[0] - dppar[2] - 2.0 / 3.0);
    coef[2] = 4.0 * cpar[0] * (cpar[0] * cpar[0] - dppar[2]);
    coef[3] = (cpar[0] * cpar[0] - dppar[2]) * (cpar[0] * cpar[0] - dppar[2])
              - dppar[0] * dppar[0];

    feff::math::ccubic(coef, sol, nsol);

    // Test solutions: must be real and within [limit1, limit2]
    for (int i = 0; i < nsol; ++i) {
        double test = std::abs(
            (cpar[0] + sol[i]) * (cpar[0] + sol[i]) - dppar[2]
            + std::sqrt(sol[i] * sol[i] * sol[i] * sol[i]
                        + 4.0 / 3.0 * sol[i] * sol[i]
                        + dppar[0] * dppar[0]));
        if (test < zero_tol) {
            if (sol[i].real() >= limit1.real() &&
                sol[i].real() <= limit2.real() &&
                std::abs(sol[i].imag()) <= zero_tol) {
                xsing[nsing] = sol[i].real();
                ++nsing;
            }
        }
    }

    // Now solve eq. 1 for q with - sign
    coef[0] = -coef[0];
    coef[2] = -coef[2];

    feff::math::ccubic(coef, sol, nsol);

    // Test solutions
    for (int i = 0; i < nsol; ++i) {
        double test = std::abs(
            (cpar[0] - sol[i]) * (cpar[0] - sol[i]) - dppar[2]
            - std::sqrt(sol[i] * sol[i] * sol[i] * sol[i]
                        + 4.0 / 3.0 * sol[i] * sol[i]
                        + dppar[0] * dppar[0]));
        if (test < zero_tol) {
            if (sol[i].real() >= limit1.real() &&
                sol[i].real() <= limit2.real() &&
                std::abs(sol[i].imag()) <= zero_tol) {
                xsing[nsing] = sol[i].real();
                ++nsing;
            }
        }
    }

    // If ifcn == 1 (singularities of r1), also solve eq. 2:
    //   q^4 + 4/3*q^2 + Wp^2 - (1 - E)^2 = 0
    if (ifcn == 1) {
        FeffComplex coef2[3];
        coef2[0] = 1.0;
        coef2[1] = 4.0 / 3.0;
        coef2[2] = dppar[0] * dppar[0];

        feff::math::cqdrtc(coef2, sol, nsol);

        for (int i = 0; i < nsol; ++i) {
            xsing2[2 * i]     =  std::sqrt(sol[i]).real();
            xsing2[2 * i + 1] = -std::sqrt(sol[i]).real();
        }

        // Test solutions
        for (int i = 0; i < 2 * nsol; ++i) {
            if (xsing2[i].real() >= limit1.real() &&
                xsing2[i].real() <= limit2.real() &&
                std::abs(sol[i].imag()) <= zero_tol) {
                xsing[nsing] = xsing2[i];
                ++nsing;
            }
        }
    }

    // Sort singularities by real part (ascending)
    if (nsing > 0) {
        std::vector<double> real_parts(nsing);
        for (int i = 0; i < nsing; ++i) {
            real_parts[i] = xsing[i].real();
        }
        auto order = feff::common::argsort(real_parts);
        std::vector<FeffComplex> temp(nsing);
        for (int i = 0; i < nsing; ++i) {
            temp[i] = xsing[i];
        }
        for (int i = 0; i < nsing; ++i) {
            xsing[i] = temp[order[i]];
        }
    }
}

} // namespace feff::exch
