// Convolution with excitation spectrum.
// Converted from: src/FF2X/exconv.f
//
// Uses the fact that for e(i+1), the convolution integral with exp(-e/ed)
// for e < e(i) is simply scaled by exp((e(i)-e(i+1))/ed). This makes
// the convolution fast.
// Written by A. Ankudinov, December 1995.

#include "exconv.hpp"
#include "ff2x_types.hpp"

#include "../math/interpolation.hpp"

#include <cmath>
#include <stdexcept>

namespace feff::ff2x {

void exconv(const double omega[], int nk, double efermi,
            double s02, double erelax, double wp, double xmu[]) {

    if (s02 >= 0.999) return;
    if (wp <= 0.0) wp = 0.00001;
    if (nk > nfinex) {
        throw std::runtime_error("check nfinex in exconv");
    }

    // Weight for plasmon excitation
    double wwp = 0.00;
    // sm1 = weight for shake up/off processes
    double sm1 = 1.0 - s02 - wwp;
    double edp = wp;
    double ed = (erelax - wwp * (wp + edp)) / sm1;

    int i0 = math::locat(efermi, nk, omega);

    double slope[nfinex]{};
    double dmu[nfinex]{};
    double xmup[nfinex]{};

    for (int i = 0; i <= i0; ++i) {
        slope[i] = 0.0;
        dmu[i] = 0.0;
    }
    for (int i = i0; i < nk - 1; ++i) {
        slope[i] = ed * (xmu[i + 1] - xmu[i]) / (omega[i + 1] - omega[i]);
    }

    double xmuf = 0.0;
    math::terp(omega, xmu, nk, 1, efermi, xmuf);

    // Start induction
    double xmult = std::exp((efermi - omega[i0 + 1]) / ed);
    dmu[i0 + 1] = xmu[i0 + 1] - slope[i0] - xmult * (xmuf - slope[i0]);
    for (int i = i0 + 1; i < nk - 1; ++i) {
        xmult = std::exp((omega[i] - omega[i + 1]) / ed);
        dmu[i + 1] = xmu[i + 1] - slope[i] + xmult * (dmu[i] - xmu[i] + slope[i]);
    }
    for (int i = 0; i < nk; ++i) {
        xmup[i] = s02 * xmu[i] + sm1 * dmu[i];
    }

    // Convolution with plasmon pole
    for (int i = i0; i < nk - 1; ++i) {
        slope[i] = slope[i] / ed * edp;
    }
    xmult = std::exp((efermi - omega[i0 + 1]) / edp);
    dmu[i0 + 1] = xmu[i0 + 1] - slope[i0] - xmult * (xmuf - slope[i0]);
    for (int i = i0 + 1; i < nk - 1; ++i) {
        xmult = std::exp((omega[i] - omega[i + 1]) / edp);
        dmu[i + 1] = xmu[i + 1] - slope[i] + xmult * (dmu[i] - xmu[i] + slope[i]);
    }

    for (int i = 0; i < nk; ++i) {
        double en = omega[i] - wp;
        int j0 = math::locat(en, nk, omega);
        if (en > efermi) {
            xmult = std::exp((omega[j0] - en) / edp);
            double dif = xmu[j0] - slope[j0];
            xmup[i] = xmup[i] + wwp * (xmult * (dmu[j0] - dif) + dif +
                                         slope[j0] * (en - omega[j0]) / edp);
        }
    }

    for (int i = 0; i < nk; ++i) {
        xmu[i] = xmup[i];
    }
}

} // namespace feff::ff2x
