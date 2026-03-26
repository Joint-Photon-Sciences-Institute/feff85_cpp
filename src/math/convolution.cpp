#include "convolution.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <cmath>
#include <algorithm>

namespace feff::math {

// Analytical convolution kernel for linear interpolation between (x1,y1) and (x2,y2)
static FeffComplex conv1(double x1, double x2, FeffComplex y1, FeffComplex y2,
                         double x0, double xloss) {
    FeffComplex result(0.0, 0.0);

    // Real part of y
    {
        double d = (x2 - x1) / 2.0;
        double a = (y2.real() - y1.real()) / 2.0;
        double b = (y2.real() + y1.real()) / 2.0;
        FeffComplex t = d / ((x1 + x2) / 2.0 - x0 - coni * xloss);
        FeffComplex dum;
        if (std::abs(t) >= 0.1) {
            dum = 2.0 * a + (b - a / t) * std::log((1.0 + t) / (1.0 - t));
        } else {
            dum = 2.0 * b * (t + t * t * t / 3.0) - 2.0 / 3.0 * a * t * t;
        }
        result = FeffComplex(dum.imag(), 0.0);
    }

    // Imaginary part of y
    {
        double d = (x2 - x1) / 2.0;
        double a = (y2.imag() - y1.imag()) / 2.0;
        double b = (y2.imag() + y1.imag()) / 2.0;
        FeffComplex t = d / ((x1 + x2) / 2.0 - x0 - coni * xloss);
        FeffComplex dum;
        if (std::abs(t) >= 0.1) {
            dum = 2.0 * a + (b - a / t) * std::log((1.0 + t) / (1.0 - t));
        } else {
            dum = 2.0 * b * (t + t * t * t / 3.0) - 2.0 / 3.0 * a * t * t;
        }
        result += coni * dum.imag();
    }

    return result;
}

void conv(const double omega[], FeffComplex xsec[], int ne1, double vicorr) {
    FeffComplex xsec0[nex];

    for (int ie = 0; ie < ne1; ++ie) {
        xsec0[ie] = FeffComplex(0.0, 0.0);
        double omega0 = omega[ie];

        // Extra point correction for finite grid at large energies
        double dx = std::max(omega[ne1 - 1] - omega[ne1 - 2], 50.0 * vicorr);
        double xlast = omega[ne1 - 1] + dx;
        dx = dx / (omega[ne1 - 1] - omega[ne1 - 2]);
        FeffComplex xsecdx = xsec[ne1 - 1] + (xsec[ne1 - 1] - xsec[ne1 - 2]) * dx;

        // First intervals
        for (int i = 0; i < ne1 - 1; ++i) {
            xsec0[ie] += conv1(omega[i], omega[i + 1], xsec[i], xsec[i + 1],
                               omega0, vicorr);
        }
        // Last interval
        xsec0[ie] += conv1(omega[ne1 - 1], xlast, xsec[ne1 - 1], xsecdx,
                           omega0, vicorr);
        xsec0[ie] /= pi;
    }

    for (int ie = 0; ie < ne1; ++ie) {
        xsec[ie] = xsec0[ie];
    }
}

} // namespace feff::math
