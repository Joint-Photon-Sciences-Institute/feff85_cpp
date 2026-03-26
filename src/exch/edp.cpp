#include "edp.hpp"
#include <feff/constants.hpp>
#include <cmath>

namespace feff::exch {

void edp(double rs, double xk, double& vr) {
    vr = 0.0;
    if (rs <= 100.0) {
        // p = sqrt(k^2 + kf^2) is the local momentum, x = p / kf
        // Reference formula 23 in "Role of Inelastic effects in EXAFS"
        // by Rehr and Chou, EXAFS1 conference edited by Bianconi.
        double xf = feff::fa / rs;
        double x = xk / xf;
        x = x + 1.0e-5;
        // Set to Fermi level if below Fermi level
        if (x < 1.00001) x = 1.00001;
        double c = std::abs((1.0 + x) / (1.0 - x));
        c = std::log(c);
        vr = -(xf / feff::pi) * (1.0 + c * (1.0 - x * x) / (2.0 * x));
    }
}

} // namespace feff::exch
