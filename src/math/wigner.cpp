#include "wigner.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace feff::math {

// Log-factorial table, shared between cwig3j and rotwig
static constexpr int idim = 58;
static double al[idim + 1];
static bool al_initialized = false;

static void init_log_factorials() {
    if (al_initialized) return;
    al[0] = 0.0;
    for (int i = 1; i <= idim; ++i) {
        al[i] = al[i - 1] + std::log(static_cast<double>(i));
    }
    al_initialized = true;
}

double cwig3j(int j1, int j2, int j3, int m1, int m2, int ient) {
    init_log_factorials();

    int m3 = -m1 - m2;
    double result = 0.0;

    if ((ient - 1) * (ient - 2) != 0) {
        throw std::runtime_error("error in cwig3j: invalid ient");
    }

    int ii = ient + ient;

    // Test parity
    if ((std::abs(m1) + std::abs(m2)) == 0 && (j1 + j2 + j3) % ii != 0) {
        return 0.0;
    }

    int m[12];
    m[0] = j1 + j2 - j3;
    m[1] = j2 + j3 - j1;
    m[2] = j3 + j1 - j2;
    m[3] = j1 + m1;
    m[4] = j1 - m1;
    m[5] = j2 + m2;
    m[6] = j2 - m2;
    m[7] = j3 + m3;
    m[8] = j3 - m3;
    m[9] = j1 + j2 + j3 + ient;
    m[10] = j2 - j3 - m1;
    m[11] = j1 - j3 + m2;

    for (int i = 0; i < 12; ++i) {
        if (i < 10) {
            if (m[i] < 0) return 0.0;
        }
        if (m[i] % ient != 0) {
            throw std::runtime_error("error in cwig3j: non-integer argument");
        }
        m[i] = m[i] / ient;
        if (i < 10 && m[i] > idim) {
            throw std::runtime_error("error in cwig3j: argument too large");
        }
    }

    int max0 = std::max({m[10], m[11], 0}) + 1;
    int min0 = std::min({m[0], m[4], m[5]}) + 1;

    int isig = 1;
    if ((max0 - 1) % 2 != 0) isig = -isig;

    double c = -al[m[9]];
    for (int i = 0; i < 9; ++i) {
        c += al[m[i]];
    }
    c /= 2.0;

    for (int i = max0; i <= min0; ++i) {
        int j = 2 - i;
        double b = al[i - 1] + al[j + m[0] - 1] + al[j + m[4] - 1] +
                   al[j + m[5] - 1] + al[i - m[10] - 1] + al[i - m[11] - 1];
        result += isig * std::exp(c - b);
        isig = -isig;
    }

    if ((j1 - j2 - m3) % ii != 0) result = -result;

    return result;
}

double rotwig(double beta, int jj, int m1, int m2, int ient) {
    init_log_factorials();

    if ((ient - 1) * (ient - 2) != 0) {
        throw std::runtime_error("Illegal ient in rotwig");
    }

    int m1p, m2p, isign;
    double betap;

    if (m1 >= 0 && std::abs(m1) >= std::abs(m2)) {
        m1p = m1;
        m2p = m2;
        betap = beta;
        isign = 1;
    } else if (m2 >= 0 && std::abs(m2) >= std::abs(m1)) {
        m1p = m2;
        m2p = m1;
        betap = -beta;
        isign = 1;
    } else if (m1 <= 0 && std::abs(m1) >= std::abs(m2)) {
        m1p = -m1;
        m2p = -m2;
        betap = beta;
        int exp = (m1 - m2) / ient;
        isign = (exp % 2 == 0) ? 1 : -1;
    } else {
        m1p = -m2;
        m2p = -m1;
        betap = -beta;
        int exp = (m2 - m1) / ient;
        isign = (exp % 2 == 0) ? 1 : -1;
    }

    double temp = 0.0;
    double zeta = std::cos(betap / 2.0);
    double eta = std::sin(betap / 2.0);

    for (int it = m1p - m2p; it <= jj - m2p; it += ient) {
        int mv[10];
        mv[0] = 1 + (jj + m1p) / ient;
        mv[1] = 1 + (jj - m1p) / ient;
        mv[2] = 1 + (jj + m2p) / ient;
        mv[3] = 1 + (jj - m2p) / ient;
        mv[4] = 1 + (jj + m1p - it) / ient;
        mv[5] = 1 + (jj - m2p - it) / ient;
        mv[6] = 1 + it / ient;
        mv[7] = 1 + (m2p - m1p + it) / ient;
        int pow_zeta = (2 * jj + m1p - m2p - 2 * it) / ient;
        int pow_eta = (2 * it - m1p + m2p) / ient;

        double factor = 0.0;
        for (int i = 0; i < 4; ++i) {
            factor += al[mv[i] - 1] / 2.0 - al[mv[i + 4] - 1];
        }

        int it_sign = ((it / ient) % 2 == 0) ? 1 : -1;

        if (pow_eta == 0 && pow_zeta == 0) {
            temp += it_sign * std::exp(factor);
        } else if (pow_eta == 0) {
            temp += it_sign * std::pow(zeta, pow_zeta) * std::exp(factor);
        } else if (pow_zeta == 0) {
            temp += it_sign * std::pow(eta, pow_eta) * std::exp(factor);
        } else {
            temp += it_sign * std::pow(zeta, pow_zeta) * std::pow(eta, pow_eta) * std::exp(factor);
        }
    }

    return isign * temp;
}

} // namespace feff::math
