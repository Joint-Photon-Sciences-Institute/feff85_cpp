#include "phase_amplitude.hpp"
#include <feff/constants.hpp>
#include <cmath>

namespace feff::math {

void atancc(FeffComplex temp, FeffComplex& phx) {
    double xx = temp.real();
    double yy = temp.imag();

    double alph;
    if (xx != 0.0) {
        double val = 1.0 - xx * xx - yy * yy;
        alph = std::sqrt(val * val + 4.0 * xx * xx) - val;
        alph = alph / (2.0 * xx);
        alph = std::atan(alph);
    } else {
        alph = 0.0;
    }

    double beta = (xx * xx + (yy + 1.0) * (yy + 1.0)) /
                  (xx * xx + (yy - 1.0) * (yy - 1.0));
    beta = std::log(beta) / 4.0;
    phx = FeffComplex(alph, beta);
}

void atan2c(FeffComplex a, FeffComplex b, FeffComplex& ampl, FeffComplex& phx) {
    double aa = std::abs(a);
    double bb = std::abs(b);

    if (aa + bb == 0.0) {
        ampl = FeffComplex(0.0, 0.0);
        phx = FeffComplex(0.0, 0.0);
    } else if (aa > bb) {
        FeffComplex temp = b / a;
        atancc(temp, phx);
        ampl = a / std::cos(phx);
    } else {
        FeffComplex temp = a / b;
        atancc(temp, phx);
        phx = pi / 2.0 - phx;
        ampl = b / std::sin(phx);
    }

    if (ampl.real() < 0.0) {
        ampl = -ampl;
        phx = phx + pi;
    }
}

void phamp(double rmt, FeffComplex pu, FeffComplex qu, FeffComplex ck,
           FeffComplex jl, FeffComplex nl, FeffComplex jlp, FeffComplex nlp,
           int ikap, FeffComplex& ph, FeffComplex& amp) {
    FeffComplex xkr = ck * rmt;
    int isign = (ikap < 0) ? -1 : 1;
    FeffComplex a_ck = ck * alphfs;
    FeffComplex factor = static_cast<double>(isign) * a_ck / (1.0 + std::sqrt(1.0 + a_ck * a_ck));

    // Find a and b such that pu = rmt*(a*jl + b*nl), qu = factor*rmt*(a*jlp + b*nlp)
    FeffComplex a = static_cast<double>(isign) * ck * xkr * (pu * nlp - qu * nl / factor);
    FeffComplex b = static_cast<double>(isign) * ck * xkr * (qu * jl / factor - pu * jlp);

    // tan(ph) = -b/a
    b = -b;
    atan2c(a, b, amp, ph);
}

} // namespace feff::math
