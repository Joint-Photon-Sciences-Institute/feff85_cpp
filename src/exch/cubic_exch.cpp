#include "cubic_exch.hpp"
#include <cmath>
#include <complex>

namespace feff::exch {

void cubic_exch(double xk0, double wp, double alph,
                double& rad, double& qplus, double& qminus) {
    constexpr double third = 1.0 / 3.0;

    // This subroutine finds the roots of the equation
    // 4*xk0 * q^3 + (alph - 4*xk0^2) * q^2 + wp^2 = 0
    // See Abramowitz and Stegun pg 17 for formulae.

    double a2 = (alph / (4.0 * xk0 * xk0) - 1.0) * xk0;
    double a0 = wp * wp / (4.0 * xk0);
    double a1 = 0.0;
    double q = a1 / 3.0 - a2 * a2 / 9.0;
    double r = (a1 * a2 - 3.0 * a0) / 6.0 - a2 * a2 * a2 / 27.0;
    rad = q * q * q + r * r;

    if (rad > 0.0) {
        qplus = 0.0;
        qminus = 0.0;
        return;
    }

    std::complex<double> s13(r, std::sqrt(-rad));
    std::complex<double> s1 = std::pow(s13, third);
    double qz1 = 2.0 * s1.real() - a2 / 3.0;
    double qz3 = -(s1.real() - std::sqrt(3.0) * s1.imag() + a2 / 3.0);
    qplus = qz1;
    qminus = qz3;
}

} // namespace feff::exch
