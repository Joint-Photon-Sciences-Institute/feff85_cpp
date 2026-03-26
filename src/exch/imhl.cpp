#include "imhl.hpp"
#include "ffq.hpp"
#include "cubic_exch.hpp"
#include "quinn.hpp"
#include <feff/constants.hpp>
#include <cmath>

namespace feff::exch {

void imhl(double rs, double xk, double& eim, int& icusp) {
    // alph is Hedin-Lundqvist parameter
    constexpr double alph = 4.0 / 3.0;

    icusp = 0;
    double xf = feff::fa / rs;
    double ef = xf * xf / 2.0;

    // xk0 is xk normalized by k Fermi
    double xk0 = xk / xf;
    // Set to Fermi level if below Fermi level
    if (xk0 < 1.00001) {
        xk0 = 1.00001;
    }

    // wp is given in units of the Fermi energy
    double wp = std::sqrt(3.0 / (rs * rs * rs)) / ef;
    double xs = wp * wp - (xk0 * xk0 - 1.0) * (xk0 * xk0 - 1.0);

    eim = 0.0;
    if (xs < 0.0) {
        double q2 = std::sqrt((std::sqrt(alph * alph - 4.0 * xs) - alph) / 2.0);
        double qu = std::min(q2, 1.0 + xk0);
        double d1 = qu - (xk0 - 1.0);
        if (d1 > 0.0) {
            eim = ffq(qu, ef, xk, wp, alph) - ffq(xk0 - 1.0, ef, xk, wp, alph);
        }
    }

    double rad, qplus, qminus;
    cubic_exch(xk0, wp, alph, rad, qplus, qminus);

    if (rad <= 0.0) {
        double d2 = qplus - (xk0 + 1.0);
        if (d2 > 0.0) {
            eim = eim + ffq(qplus, ef, xk, wp, alph) -
                        ffq(xk0 + 1.0, ef, xk, wp, alph);
        }
        double d3 = (xk0 - 1.0) - qminus;
        if (d3 > 0.0) {
            eim = eim + ffq(xk0 - 1.0, ef, xk, wp, alph) -
                        ffq(qminus, ef, xk, wp, alph);
            // Beginning of the imaginary part and position of the cusp x0
            icusp = 1;
        }
    }

    double ei;
    quinn(xk0, rs, wp, ef, ei);
    if (eim >= ei) eim = ei;
}

} // namespace feff::exch
