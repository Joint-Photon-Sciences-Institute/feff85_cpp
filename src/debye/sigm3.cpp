// Third cumulant from Einstein model with Morse potential.
// Converted from: src/DEBYE/sigm3.f
// Reference: Nguyen Van Hung & J.J.Rehr, Phys. Rev. B 56, 43

#include "debye.hpp"

#include <feff/constants.hpp>

#include <cmath>

namespace feff::debye {

void sigm3(double& sig1, double sig2, double& sig3,
           double tk, double alphat, double thetae) {

    constexpr double four_thirds = 4.0 / 3.0;
    constexpr double three_quarters = 3.0 / 4.0;

    // Convert alphat to code units (multiply by bohr)
    double alphat_cu = alphat * bohr;

    double z = std::exp(-thetae / tk);
    double sig02 = (1.0 - z) / (1.0 + z) * sig2;
    double sig01 = alphat_cu * sig02 * three_quarters;
    sig1 = sig01 * sig2 / sig02;
    sig3 = (2.0 - four_thirds * (sig02 / sig2) * (sig02 / sig2)) * sig1 * sig2;
}

} // namespace feff::debye
