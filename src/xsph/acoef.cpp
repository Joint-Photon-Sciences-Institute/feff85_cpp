// Angular coefficient matrix for spin/orbital moment calculations.
// Converted from src/XSPH/acoef.f

#include "acoef.hpp"
#include <feff/dimensions.hpp>
#include <cmath>
#include <cstdio>

namespace feff::math { double cwig3j(int j1, int j2, int j3, int m1, int m2, int m3); }

namespace feff::xsph {

void kfromi(int i, int lpp, int& jj, int& k) {
    jj = lpp + i - 2;
    k = -lpp - 1;
    if (i == 1) k = lpp;
}

void acoef(int ispin, float amat[]) {
    // Dimensions: (2*lx+1) * 2 * 2 * 3 * (lx+1)
    int total = (2 * lx + 1) * 2 * 2 * 3 * (lx + 1);
    for (int i = 0; i < total; i++) amat[i] = 0.0f;

    int ms = 1;
    if (ispin < 0) ms = 0;
    std::printf(" Spin = %d\n", 2 * ms - 1);

    // Tabulate Clebsch-Gordan coefficients
    // t3j[lp][jp][mp] dimensions: (lx+1) x (lx+1) x 2
    float t3j[(lx + 1)][(lx + 1)][2];
    float operls[2][2][3];

    for (int ml = -lx; ml <= lx; ml++) {
        int mj = 2 * ml + (2 * ms - 1);
        float xmj = 0.5f * mj;
        mj = -mj;

        // Tabulate necessary CG coefficients
        for (int lp = 0; lp <= lx; lp++) {
            for (int jp = 0; jp <= lx; jp++) {
                for (int mp = 0; mp < 2; mp++) {
                    int lp2 = 2 * lp;
                    int jp2 = 2 * jp + 1;
                    int mp2 = 2 * mp - 1;
                    int sign = (lp % 2 == 0) ? 1 : -1;
                    t3j[lp][jp][mp] = sign * std::sqrt(static_cast<float>(jp2 + 1)) *
                        static_cast<float>(feff::math::cwig3j(1, jp2, lp2, mp2, mj, 2));
                }
            }
        }

        for (int lpp = 0; lpp <= lx; lpp++) {
            for (int m1 = 0; m1 < 2; m1++) {
                for (int m2 = 0; m2 < 2; m2++) {
                    for (int iop = 0; iop < 3; iop++) {
                        operls[m1][m2][iop] = 0.0f;
                        if (m1 == m2) {
                            float xms = m1 - 0.5f;
                            float xml = xmj - xms;
                            if (std::abs(ml + ms - m1) <= lpp) {
                                if (ispin == 0) {
                                    operls[m1][m2][iop] = 2.0f;
                                } else if (iop == 0) {
                                    operls[m1][m2][iop] = xms;
                                } else if (iop == 1 && std::abs(ispin) == 1) {
                                    operls[m1][m2][iop] = xml;
                                } else if (iop == 1 && std::abs(ispin) == 2) {
                                    operls[m1][m2][iop] = 1.0f;
                                } else if (iop == 2 && std::abs(ispin) == 1) {
                                    operls[m1][m2][iop] = xms * 2.0f *
                                        (3.0f * xml * xml - lpp * (lpp + 1)) /
                                        (2 * lpp + 3) / (2 * lpp - 1);
                                } else if (iop == 2 && std::abs(ispin) == 2) {
                                    operls[m1][m2][iop] = t3j[lpp][lpp][m1] * t3j[lpp][lpp][m1];
                                }
                            }
                        } else {
                            if (iop == 2 && std::abs(ispin) <= 1 &&
                                static_cast<int>(0.5f + std::abs(xmj)) < lpp) {
                                operls[m1][m2][iop] = 3.0f * xmj *
                                    std::sqrt(lpp * (lpp + 1) - (xmj * xmj - 0.25f)) /
                                    (2 * lpp + 3) / (2 * lpp - 1);
                            } else if (iop == 2 && std::abs(ispin) > 1) {
                                operls[m1][m2][iop] = t3j[lpp][lpp][m1] * t3j[lpp][lpp][m2];
                            }
                        }
                    }
                }
            }

            // Calculate amat
            for (int i1 = 0; i1 < 2; i1++) {
                int jj, k1;
                kfromi(i1 + 1, lpp, jj, k1);
                if (k1 != 0) {
                    for (int i2 = 0; i2 < 2; i2++) {
                        int jp, k2;
                        kfromi(i2 + 1, lpp, jp, k2);
                        if (k2 != 0) {
                            for (int iop = 0; iop < 3; iop++) {
                                for (int m2 = 0; m2 < 2; m2++) {
                                    for (int m1 = 0; m1 < 2; m1++) {
                                        amat[amat_index(ml, i1, i2, iop, lpp)] +=
                                            operls[m1][m2][iop] *
                                            t3j[lpp][jp][ms] * t3j[lpp][jp][m1] *
                                            t3j[lpp][jj][m2] * t3j[lpp][jj][ms];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

} // namespace feff::xsph
