// Third cumulant from Debye model (single scattering only).
// Converted from: src/DEBYE/sigte3.f
// Uses reduced mass and Einstein frequency to compute sig1, sig3.

#include "debye.hpp"

#include "../common/periodic_table.hpp"

#include <cmath>

namespace feff::debye {

// Physical constants (updated 2017)
static constexpr double hbar_si = 1.054571800e-34;  // J*s
static constexpr double amu_si  = 1.660539040e-27;  // kg
static constexpr double xkb_si  = 1.38064852e-23;   // J/K

void sigte3(int iz1, int iz2, double sig2, double alphat,
            double thetad, double reff, double& sig1, double& sig3) {

    double ami = common::atomic_weights[iz1] * amu_si;
    double amj = common::atomic_weights[iz2] * amu_si;

    // Reduced mass
    double xmu = 1.0 / (1.0 / ami + 1.0 / amj);

    // Einstein frequency from Debye temperature
    double omega = (2.0 * xkb_si * thetad) / (3.0 * hbar_si);

    // Spring constants
    double xks = xmu * omega * omega;
    double xk3 = xks * xks * static_cast<double>(reff) * alphat / (3.0 * xkb_si);

    // Zero-point sigma^2
    double sig02 = hbar_si * omega / xks;

    // First cumulant
    sig1 = -3.0 * (xk3 / xks) * sig2;

    // Third cumulant
    sig3 = 2.0 - (4.0 / 3.0) * (sig02 / sig2) * (sig02 / sig2);
    sig3 = sig3 * (sig1 * sig2);
}

} // namespace feff::debye
