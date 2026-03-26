#include "bcoef.hpp"
#include "wigner.hpp"
#include <cmath>
#include <cstring>

namespace feff::math {

void bcoef(int kinit, int ipol, const FeffComplex ptz[3][3], int le2,
           bool ltrace, int ispin, double angks,
           int kind[8], int lind[8],
           FeffComplex bmat_flat[(2*lx+1)*2*8*(2*lx+1)*2*8]) {

    // Helper lambda for 6D bmat access: bmat(ml1, ms1, k1, ml2, ms2, k2)
    auto bmat = [&](int ml1, int ms1, int k1, int ml2, int ms2, int k2) -> FeffComplex& {
        return bmat_flat[bmat_index(ml1, ms1, k1, ml2, ms2, k2)];
    };

    // Zero bmat
    int total_size = (2 * lx + 1) * 2 * 8 * (2 * lx + 1) * 2 * 8;
    for (int i = 0; i < total_size; ++i) {
        bmat_flat[i] = FeffComplex(0.0, 0.0);
    }

    int jind[8];

    // 3 dipole transitions
    for (int k = -1; k <= 1; ++k) {
        int kap = kinit + k;
        if (k == 0) kap = -kap;
        int jkap = std::abs(kap);
        int lkap = kap;
        if (kap <= 0) lkap = std::abs(kap) - 1;
        if (lkap > lx) {
            jkap = 0;
            lkap = -1;
            kap = 0;
        }
        jind[k + 1] = jkap;
        lind[k + 1] = lkap;
        kind[k + 1] = kap;
    }

    // 5 quadrupole or 3 mag.dipole transitions
    for (int k = -2; k <= 2; ++k) {
        int jkap = std::abs(kinit) + k;
        if (jkap <= 0) jkap = 0;
        int kap = jkap;
        if (kinit < 0 && std::abs(k) != 1) kap = -jkap;
        if (kinit > 0 && std::abs(k) == 1) kap = -jkap;
        int lkap = kap;
        if (kap <= 0) lkap = -kap - 1;
        if (lkap > lx || le2 == 0 || (le2 == 1 && std::abs(k) == 2)) {
            jkap = 0;
            lkap = -1;
            kap = 0;
        }
        jind[k + 5] = jkap;
        lind[k + 5] = lkap;
        kind[k + 5] = kap;
    }

    if (ipol == 0) {
        // Polarization average: bmat is diagonal and simple
        for (int k = 0; k < 8; ++k) {
            for (int ms = 0; ms <= 1; ++ms) {
                for (int ml = -lind[k]; ml <= lind[k]; ++ml) {
                    int i2 = 3;
                    if (le2 == 2 && k > 2) i2 = 5;
                    FeffComplex val(0.5 / (2.0 * lind[k] + 1.0) / i2, 0.0);
                    if (k < 3) val = -val;
                    bmat(ml, ms, k, ml, ms, k) = val;
                }
            }
        }
    } else {
        // Full polarization calculation (ipol=1 or ipol=2)

        // t3j and x3j arrays
        // t3j[k][ms][mp], x3j[k][p][mp] where mp ranges [-lx, lx+1]
        constexpr int mp_dim = 2 * lx + 2;  // -lx to lx+1
        auto mp_idx = [](int mp) { return mp + lx; };

        double t3j[8][2][mp_dim];
        double x3j[8][3][mp_dim];  // p index: -1,0,1 → 0,1,2
        std::memset(t3j, 0, sizeof(t3j));
        std::memset(x3j, 0, sizeof(x3j));

        for (int k1 = 0; k1 < 8; ++k1) {
            for (int mp = -jind[k1] + 1; mp <= jind[k1]; ++mp) {
                for (int ms = 0; ms <= 1; ++ms) {
                    int j1 = 2 * lind[k1];
                    int j2 = 1;
                    int j3 = 2 * jind[k1] - 1;
                    int m1 = 2 * (mp - ms);
                    int m2 = 2 * ms - 1;
                    t3j[k1][ms][mp_idx(mp)] = std::sqrt(j3 + 1.0) *
                        cwig3j(j1, j2, j3, m1, m2, 2);
                    if ((j2 - j1 - m1 - m2) / 2 % 2 != 0)
                        t3j[k1][ms][mp_idx(mp)] = -t3j[k1][ms][mp_idx(mp)];
                }
                for (int i = -1; i <= 1; ++i) {
                    int j1 = 2 * jind[k1] - 1;
                    int j2 = 2;
                    if (k1 > 2 && le2 == 2) j2 = 4;
                    int j3 = 2 * std::abs(kinit) - 1;
                    int m1 = -2 * mp + 1;
                    int m2 = 2 * i;
                    x3j[k1][i + 1][mp_idx(mp)] = cwig3j(j1, j2, j3, m1, m2, 2);
                }
            }
        }

        // qmat[mj][ml][ms][k]
        constexpr int mj_dim = mp_dim;
        double qmat[mj_dim][(2*lx+1)][2][8];
        std::memset(qmat, 0, sizeof(qmat));

        for (int i = 0; i < 8; ++i) {
            for (int ms = 0; ms <= 1; ++ms) {
                for (int ml = -lind[i]; ml <= lind[i]; ++ml) {
                    for (int mj = -jind[i] + 1; mj <= jind[i]; ++mj) {
                        int mp = ml + ms;
                        int jj = 2 * jind[i] - 1;
                        int mmj = 2 * mj - 1;
                        int mmp = 2 * mp - 1;
                        double value = rotwig(angks, jj, mmj, mmp, 2);
                        qmat[mp_idx(mj)][ml + lx][ms][i] =
                            value * t3j[i][ms][mp_idx(mp)];
                    }
                }
            }
        }

        // pmat[m1][i1][m2][i2]
        FeffComplex pmat[mj_dim][8][mj_dim][8];
        std::memset(pmat, 0, sizeof(pmat));

        for (int i2 = 0; i2 < 8; ++i2) {
            for (int m2 = -jind[i2] + 1; m2 <= jind[i2]; ++m2) {
                for (int i1 = 0; i1 < 8; ++i1) {
                    for (int m1 = -jind[i1] + 1; m1 <= jind[i1]; ++m1) {
                        if (std::abs(m2 - m1) <= 2) {
                            FeffComplex sum(0.0, 0.0);
                            for (int j = -1; j <= 1; ++j) {
                                for (int ii = -1; ii <= 1; ++ii) {
                                    if (m1 - ii == m2 - j) {
                                        int is = 1;
                                        if (le2 == 1 && ii > 0 && i1 > 2) is = -is;
                                        if (le2 == 1 && j > 0 && i2 > 2) is = -is;
                                        sum += static_cast<double>(is) *
                                               x3j[i1][ii + 1][mp_idx(m1)] *
                                               ptz[ii + 1][j + 1] *
                                               x3j[i2][j + 1][mp_idx(m2)];
                                    }
                                }
                            }
                            // Phase factor
                            int is = 1;
                            if ((jind[i1] - jind[i2]) % 2 != 0) is = -is;
                            if (i2 < 3) is = -is;
                            FeffComplex phase = std::pow(coni, lind[i2] - lind[i1]);
                            pmat[mp_idx(m1)][i1][mp_idx(m2)][i2] =
                                sum * static_cast<double>(is) * phase;
                        }
                    }
                }
            }
        }

        // tmat = pmat * qmat
        FeffComplex tmat[mj_dim][8][(2*lx+1)][2][8];
        std::memset(tmat, 0, sizeof(tmat));

        for (int i1 = 0; i1 < 8; ++i1) {
            for (int ms = 0; ms <= 1; ++ms) {
                for (int ml = -lind[i1]; ml <= lind[i1]; ++ml) {
                    for (int i2 = 0; i2 < 8; ++i2) {
                        for (int mj = -jind[i2] + 1; mj <= jind[i2]; ++mj) {
                            FeffComplex sum(0.0, 0.0);
                            for (int mp = -jind[i1] + 1; mp <= jind[i1]; ++mp) {
                                sum += pmat[mp_idx(mj)][i2][mp_idx(mp)][i1] *
                                       qmat[mp_idx(mp)][ml + lx][ms][i1];
                            }
                            tmat[mp_idx(mj)][i2][ml + lx][ms][i1] = sum;
                        }
                    }
                }
            }
        }

        // bmat = qmat^T * tmat
        for (int i1 = 0; i1 < 8; ++i1) {
            for (int ms1 = 0; ms1 <= 1; ++ms1) {
                for (int ml1 = -lind[i1]; ml1 <= lind[i1]; ++ml1) {
                    for (int i2 = 0; i2 < 8; ++i2) {
                        for (int ms2 = 0; ms2 <= 1; ++ms2) {
                            for (int ml2 = -lind[i2]; ml2 <= lind[i2]; ++ml2) {
                                FeffComplex sum(0.0, 0.0);
                                for (int mj = -jind[i2] + 1; mj <= jind[i2]; ++mj) {
                                    sum += qmat[mp_idx(mj)][ml2 + lx][ms2][i2] *
                                           tmat[mp_idx(mj)][i2][ml1 + lx][ms1][i1];
                                }
                                bmat(ml2, ms2, i2, ml1, ms1, i1) = sum;
                            }
                        }
                    }
                }
            }
        }
    }

    // Trace over ml if needed
    if (ltrace) {
        for (int i1 = 0; i1 < 8; ++i1) {
            for (int ms1 = 0; ms1 <= 1; ++ms1) {
                for (int i2 = 0; i2 < 8; ++i2) {
                    for (int ms2 = 0; ms2 <= 1; ++ms2) {
                        if (lind[i1] != lind[i2] || ms1 != ms2) {
                            bmat(0, ms2, i2, 0, ms1, i1) = FeffComplex(0.0, 0.0);
                        } else {
                            for (int ml = 1; ml <= lind[i1]; ++ml) {
                                bmat(0, ms1, i2, 0, ms1, i1) +=
                                    bmat(-ml, ms1, i2, -ml, ms1, i1) +
                                    bmat(ml, ms1, i2, ml, ms1, i1);
                            }
                        }
                    }
                }
            }
        }
    }

    // Spin averaging / selection
    if (ispin == 0) {
        for (int i1 = 0; i1 < 8; ++i1) {
            for (int i2 = 0; i2 < 8; ++i2) {
                for (int ml1 = -lind[i1]; ml1 <= lind[i1]; ++ml1) {
                    for (int ml2 = -lind[i2]; ml2 <= lind[i2]; ++ml2) {
                        bmat(ml2, 0, i2, ml1, 0, i1) +=
                            bmat(ml2, 1, i2, ml1, 1, i1);
                    }
                }
            }
        }
    } else if (ispin == 2 || (ispin == 1 && nspx == 1)) {
        for (int i1 = 0; i1 < 8; ++i1) {
            for (int i2 = 0; i2 < 8; ++i2) {
                for (int ml1 = -lind[i1]; ml1 <= lind[i1]; ++ml1) {
                    for (int ml2 = -lind[i2]; ml2 <= lind[i2]; ++ml2) {
                        bmat(ml2, 0, i2, ml1, 0, i1) =
                            bmat(ml2, 1, i2, ml1, 1, i1);
                    }
                }
            }
        }
    }
}

} // namespace feff::math
