// Dichroism / n* calculation.
// Converted from GENFMT/xstar.f

#include "xstar.hpp"
#include <cmath>

namespace feff::genfmt {

double xxcos(const double veca[3], const double vecb[3]) {
    double x1 = 0.0;
    for (int j = 0; j < 3; ++j) {
        x1 += veca[j] * vecb[j];
    }
    double xnorma = 0.0;
    double xnormb = 0.0;
    for (int j = 0; j < 3; ++j) {
        xnorma += veca[j] * veca[j];
        xnormb += vecb[j] * vecb[j];
    }
    return x1 / std::sqrt(xnorma * xnormb);
}

double ystar(int lfin, double x, double y, double z, int iav) {
    // Legendre polynomial coefficients pln(i, lfin)
    // pln[degree][lfin-1]
    static const double pln[5][4] = {
        {  0.0,   -0.5,    0.0,    0.375 },
        {  1.0,    0.0,   -1.5,    0.0   },
        {  0.0,    1.5,    0.0,   -3.75  },
        {  0.0,    0.0,    2.5,    0.0   },
        {  0.0,    0.0,    0.0,    4.375 }
    };

    double pln0 = pln[0][lfin - 1];
    for (int i = 1; i <= lfin; ++i) {
        pln0 += pln[i][lfin - 1] * std::pow(x, i);
    }

    if (iav == 0) {
        return pln0 / (2.0 * lfin + 1.0);
    } else {
        double pln1 = pln[1][lfin - 1];
        for (int i = 2; i <= lfin; ++i) {
            pln1 += pln[i][lfin - 1] * i * std::pow(x, i - 1);
        }

        double pln2 = 2.0 * pln[2][lfin - 1];
        for (int i = 3; i <= lfin; ++i) {
            pln2 += pln[i][lfin - 1] * i * (i - 1) * std::pow(x, i - 2);
        }

        double ytemp = -lfin * pln0 + pln1 * (x + y * z) -
                        pln2 * (y * y + z * z - 2.0 * x * y * z);
        return ytemp * 3.0 / lfin / (4.0 * lfin * lfin - 1.0);
    }
}

double xstar(const double eps1[3], const double eps2[3],
             const double vec1[3], const double vec2[3],
             int ndeg, double elpty, int ilinit) {

    int lfin = ilinit;
    double x = xxcos(vec1, vec2);
    int iav = 1;
    double y = xxcos(eps1, vec1);
    double z = xxcos(eps1, vec2);
    double xtemp = ystar(lfin, x, y, z, iav);

    if (elpty != 0.0) {
        y = xxcos(eps2, vec1);
        z = xxcos(eps2, vec2);
        xtemp += elpty * elpty * ystar(lfin, x, y, z, iav);
    }

    return ndeg * xtemp / (1.0 + elpty * elpty);
}

} // namespace feff::genfmt
