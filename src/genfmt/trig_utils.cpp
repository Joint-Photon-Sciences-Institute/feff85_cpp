// Trigonometric utility functions for path geometry.
// Converted from GENFMT/trig.f

#include "trig_utils.hpp"
#include <cmath>

namespace feff::genfmt {

void trig(double x, double y, double z,
          double& ct, double& st, double& cp, double& sp) {
    constexpr double eps = 1.0e-6;

    double r = std::sqrt(x * x + y * y + z * z);
    double rxy = std::sqrt(x * x + y * y);

    if (r < eps) {
        ct = 1.0;
        st = 0.0;
    } else {
        ct = z / r;
        st = rxy / r;
    }

    if (rxy < eps) {
        cp = 1.0;
        if (ct < 0.0) cp = -1.0;
        sp = 0.0;
    } else {
        cp = x / rxy;
        sp = y / rxy;
    }
}

void arg(FeffComplex c, double fi, double& th) {
    constexpr double eps = 1.0e-6;

    double x = c.real();
    double y = c.imag();
    if (std::abs(x) < eps) x = 0.0;
    if (std::abs(y) < eps) y = 0.0;

    if (std::abs(x) < eps && std::abs(y) < eps) {
        th = fi;
    } else {
        th = std::atan2(y, x);
    }
}

} // namespace feff::genfmt
