#include "polynomial_roots.hpp"
#include <feff/constants.hpp>
#include <cmath>

namespace feff::math {

void cqdrtc(const FeffComplex coef[3], FeffComplex sol[2], int& nsol) {
    if (coef[0] == 0.0) {
        if (coef[1] == 0.0) {
            nsol = -1;
            return;
        } else {
            nsol = 1;
            sol[0] = -coef[2] / coef[1];
        }
    } else {
        nsol = 2;
        FeffComplex disc = coef[1] * coef[1] - 4.0 * coef[0] * coef[2];
        FeffComplex root = std::sqrt(disc);
        double sgn = (std::conj(coef[1]) * root).real() >= 0.0 ? 1.0 : -1.0;
        FeffComplex q = -0.5 * (coef[1] + sgn * root);
        sol[0] = q / coef[0];
        sol[1] = coef[2] / q;
    }
}

void ccubic(const FeffComplex coef[4], FeffComplex sol[4], int& nsol) {
    if (coef[0] == 0.0) {
        FeffComplex coef2[3] = {coef[1], coef[2], coef[3]};
        cqdrtc(coef2, sol, nsol);
    } else {
        FeffComplex a = coef[1] / coef[0];
        FeffComplex b = coef[2] / coef[0];
        FeffComplex c = coef[3] / coef[0];
        nsol = 3;

        FeffComplex Q = (a * a - 3.0 * b) / 9.0;
        FeffComplex R = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 54.0;

        if (Q.imag() == 0.0 && R.imag() == 0.0 &&
            (R * R).imag() < (Q * Q * Q).imag()) {
            double theta = std::acos((R / std::sqrt(Q * Q * Q)).real());
            sol[0] = -2.0 * std::sqrt(Q) * std::cos(theta / 3.0) - a / 3.0;
            sol[1] = -2.0 * std::sqrt(Q) * std::cos((theta + 2.0 * pi) / 3.0) - a / 3.0;
            sol[2] = -2.0 * std::sqrt(Q) * std::cos((theta - 2.0 * pi) / 3.0) - a / 3.0;
        } else {
            FeffComplex disc = R * R - Q * Q * Q;
            double sgn = (std::conj(R) * std::sqrt(disc)).real() >= 0.0 ? 1.0 : -1.0;
            FeffComplex P1 = -std::pow(R + sgn * std::sqrt(disc), 1.0 / 3.0);
            FeffComplex P2;
            if (P1 == 0.0) {
                P2 = FeffComplex(0.0, 0.0);
            } else {
                P2 = Q / P1;
            }
            FeffComplex I(0.0, 1.0);
            sol[0] = (P1 + P2) - a / 3.0;
            sol[1] = -0.5 * (P1 + P2) - a / 3.0 + I * std::sqrt(3.0) / 2.0 * (P1 - P2);
            sol[2] = -0.5 * (P1 + P2) - a / 3.0 - I * std::sqrt(3.0) / 2.0 * (P1 - P2);
        }
    }
}

void quartic(FeffComplex q[4]) {
    constexpr double One = 1.0, Two = 2.0, Four = 4.0, Eight = 8.0;
    constexpr double Twelve = 12.0, Twnt7 = 27.0, Sevnt2 = 72.0;
    constexpr double Tto1O3 = 1.259921049894873;   // 2^(1/3)
    constexpr double Tto2O3 = 1.587401051968199;   // 2^(2/3)
    constexpr double Root6 = 2.449489742783178;     // sqrt(6)

    FeffComplex a = q[0], b = q[1], c = q[2], d = q[3];

    FeffComplex F = b * b + Twelve * a * d;
    FeffComplex G = Two * b * b * b + Twnt7 * a * c * c - Sevnt2 * a * b * d;

    FeffComplex A1 = std::pow(G + std::sqrt(-Four * F * F * F + G * G), 1.0 / 3.0);
    FeffComplex B1 = Two * Tto1O3 * F;

    FeffComplex P = std::sqrt((-Four * b + B1 / A1 + Tto2O3 * A1) / a);

    FeffComplex D1 = Eight * b + B1 / A1 + Tto2O3 * A1;
    FeffComplex D2 = Twelve * Root6 * c / P;

    FeffComplex QPlus = std::sqrt(-(D1 + D2) / a);
    FeffComplex QMin = std::sqrt(-(D1 - D2) / a);

    double amp = One / (Two * Root6);

    q[0] = amp * (P - QPlus);
    q[1] = amp * (P + QPlus);
    q[2] = -amp * (P + QMin);
    q[3] = amp * (-P + QMin);
}

} // namespace feff::math
