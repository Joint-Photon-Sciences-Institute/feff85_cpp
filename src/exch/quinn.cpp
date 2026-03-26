#include "quinn.hpp"
#include <feff/constants.hpp>
#include <cmath>

namespace feff::exch {

void quinn(double x, double rs, double wp, double ef, double& ei) {
    // alphaq = 1/fa
    double alphaq = 1.0 / feff::fa;

    // Calculate Quinn prefactor in atomic Hartree units
    double pisqrt = std::sqrt(feff::pi);
    double pfq = pisqrt / (32.0 * std::pow(alphaq * rs, 1.5));
    double temp1 = std::atan(std::sqrt(feff::pi / (alphaq * rs)));
    double temp2 = std::sqrt(alphaq * rs / feff::pi) / (1.0 + alphaq * rs / feff::pi);
    pfq = pfq * (temp1 + temp2);

    // Calculate Quinn cutoff
    // wkc = Quinn's plasmon threshold
    // Quinn, PR126, 1453, 1962, eq. (11)
    // In formulae below wp = omegap / ef
    double wkc = (std::sqrt(1.0 + wp) - 1.0);
    wkc = wkc * wkc;
    wkc = (1.0 + (6.0 / 5.0) * wkc / (wp * wp)) * wp * ef;

    // Add Fermi energy to get correct energy for plasma excitations to turn on
    double ekc = wkc + ef;

    // Calculate gamma
    double gam = (pfq / x) * (x * x - 1.0) * (x * x - 1.0);

    // Put in Fermi function cutoff
    double eabs = ef * x * x;
    double arg = (eabs - ekc) / (0.3 * ekc);
    double f = 0.0;
    if (arg < 80.0) f = 1.0 / (1.0 + std::exp(arg));

    ei = -gam * f / 2.0;
}

} // namespace feff::exch
