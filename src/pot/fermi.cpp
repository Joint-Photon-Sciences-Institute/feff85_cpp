// Fermi level from interstitial density.
// Converted from src/POT/fermi.f

#include "fermi.hpp"
#include <feff/constants.hpp>
#include <cmath>

namespace feff::pot {

void fermi(double rhoint, double vint, double& xmu, double& rs, double& xf)
{
    // den is the interstitial density
    double den = rhoint / (4.0 * pi);

    // rs is the density parameter
    // Use cbrt() instead of pow(x, 1/3) to correctly handle any sign
    // (Fortran's ** operator with real exponent also uses real cube root)
    rs = std::cbrt(3.0 / (4.0 * pi * den));

    // xf is the interstitial Fermi momentum (kf = fa/rs)
    xf = fa / rs;

    // xmu is the Fermi level in Hartrees
    xmu = vint + xf * xf / 2.0;
}

} // namespace feff::pot
