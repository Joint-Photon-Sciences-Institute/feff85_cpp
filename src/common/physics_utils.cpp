// Physics utility functions.
// Converted from: src/COMMON/getxk.f, pijump.f, setkap.f, setgam.f, iniptz.f

#include "physics_utils.hpp"
#include "logging.hpp"
#include "../par/parallel.hpp"
#include <feff/constants.hpp>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>

// Forward declaration — terp from math module
namespace feff::math {
    void terp(const double x[], const double y[], int n, int m,
              double x0, double& y0);
}

namespace feff::common {

double getxk(double e) {
    double xk = std::sqrt(std::abs(2.0 * e));
    if (e < 0.0) xk = -xk;
    return xk;
}

void pijump(double& ph, double old) {
    constexpr double twopi = 2.0 * pi;
    double xph[3];

    xph[0] = ph - old;
    int jump = static_cast<int>((std::abs(xph[0]) + pi) / twopi);
    xph[1] = xph[0] - jump * twopi;
    xph[2] = xph[0] + jump * twopi;

    double abs0 = std::abs(xph[0]);
    double abs1 = std::abs(xph[1]);
    double abs2 = std::abs(xph[2]);
    double xphmin = std::min(abs0, std::min(abs1, abs2));

    int isave = -1;
    for (int i = 0; i < 3; ++i) {
        if (std::abs(xphmin - std::abs(xph[i])) <= 0.01) {
            isave = i;
        }
    }

    if (isave < 0) {
        feff::par::par_stop("pijump");
    }

    ph = old + xph[isave];
}

void setkap(int ihole, int& kinit, int& linit) {
    // s states: K(1), LI(2), MI(5), NI(10), OI(17), PI(24), QI(27)
    if (ihole <= 2 || ihole == 5 || ihole == 10 ||
        ihole == 17 || ihole == 24 || ihole == 27) {
        linit = 0;
        kinit = -1;
    }
    // p 1/2 states: LII(3), MII(6), NII(11), OII(18), PII(25), QII(30)
    else if (ihole == 3 || ihole == 6 || ihole == 11 ||
             ihole == 18 || ihole == 25 || ihole == 30) {
        linit = 1;
        kinit = 1;
    }
    // p 3/2 states: LIII(4), MIII(7), NIII(12), OIII(19), PIII(26)
    else if (ihole == 4 || ihole == 7 || ihole == 12 ||
             ihole == 19 || ihole == 26) {
        linit = 1;
        kinit = -2;
    }
    // d 3/2 states: MIV(8), NIV(13), OIV(20), PIV(27)
    else if (ihole == 8 || ihole == 13 ||
             ihole == 20 || ihole == 28) {
        linit = 2;
        kinit = 2;
    }
    // d 5/2 states: MV(9), NV(14), OV(21), PV(28)
    else if (ihole == 9 || ihole == 14 ||
             ihole == 21 || ihole == 29) {
        linit = 2;
        kinit = -3;
    }
    // f 5/2 states: NVI(15), OVI(22)
    else if (ihole == 15 || ihole == 22) {
        linit = 3;
        kinit = 3;
    }
    // f 7/2 states: NVII(16), OVII(23)
    else if (ihole == 16 || ihole == 23) {
        linit = 3;
        kinit = -4;
    }
    else {
        feff::par::par_stop("invalid hole number in setkap");
    }
}

void setgam(int iz, int ihole, double& gamach) {
    // Core-hole lifetime widths from:
    // K. Rahkonen and K. Krause, Atomic Data and Nuclear Data Tables, Vol 14, No 2, 1974.
    // zh[i][j] = atomic number breakpoints for hole type j
    // gamh[i][j] = gamma values (eV) at those Z values

    // Data tables: 8 Z-breakpoints x 16 hole types
    static const double zh[16][8] = {
        { 0.99, 10.0, 20.0, 40.0, 50.0, 60.0, 80.0, 95.1},  // K
        { 0.99, 18.0, 22.0, 35.0, 50.0, 52.0, 75.0, 95.1},  // LI
        { 0.99, 17.0, 28.0, 31.0, 45.0, 60.0, 80.0, 95.1},  // LII
        { 0.99, 17.0, 28.0, 31.0, 45.0, 60.0, 80.0, 95.1},  // LIII
        { 0.99, 20.0, 28.0, 30.0, 36.0, 53.0, 80.0, 95.1},  // MI
        { 0.99, 20.0, 22.0, 30.0, 40.0, 68.0, 80.0, 95.1},  // MII
        { 0.99, 20.0, 22.0, 30.0, 40.0, 68.0, 80.0, 95.1},  // MIII
        { 0.99, 36.0, 40.0, 48.0, 58.0, 76.0, 79.0, 95.1},  // MIV
        { 0.99, 36.0, 40.0, 48.0, 58.0, 76.0, 79.0, 95.1},  // MV
        { 0.99, 30.0, 40.0, 47.0, 50.0, 63.0, 80.0, 95.1},  // NI
        { 0.99, 40.0, 42.0, 49.0, 54.0, 70.0, 87.0, 95.1},  // NII
        { 0.99, 40.0, 42.0, 49.0, 54.0, 70.0, 87.0, 95.1},  // NIII
        { 0.99, 40.0, 50.0, 55.0, 60.0, 70.0, 81.0, 95.1},  // NIV
        { 0.99, 40.0, 50.0, 55.0, 60.0, 70.0, 81.0, 95.1},  // NV
        { 0.99, 71.0, 73.0, 79.0, 86.0, 90.0, 95.0, 100.0}, // NVI
        { 0.99, 71.0, 73.0, 79.0, 86.0, 90.0, 95.0, 100.0}  // NVII
    };

    static const double gamh[16][8] = {
        { 0.02,   0.28,  0.75,   4.8,   10.5,  21.0,  60.0,  105.0},
        { 0.07,   3.9,   3.8,    7.0,    6.0,   3.7,   8.0,   19.0},
        { 0.001,  0.12,  1.4,    0.8,    2.6,   4.1,   6.3,   10.5},
        { 0.001,  0.12,  0.55,   0.7,    2.1,   3.5,   5.4,    9.0},
        { 0.001,  1.0,   2.9,    2.2,    5.5,  10.0,  22.0,   22.0},
        { 0.001,  0.001, 0.5,    2.0,    2.6,  11.0,  15.0,   16.0},
        { 0.001,  0.001, 0.5,    2.0,    2.6,  11.0,  10.0,   10.0},
        { 0.0006, 0.09,  0.07,   0.48,   1.0,   4.0,   2.7,    4.7},
        { 0.0006, 0.09,  0.07,   0.48,   0.87,  2.2,   2.5,    4.3},
        { 0.001,  0.001, 6.2,    7.0,    3.2,  12.0,  16.0,   13.0},
        { 0.001,  0.001, 1.9,   16.0,    2.7,  13.0,  13.0,    8.0},
        { 0.001,  0.001, 1.9,   16.0,    2.7,  13.0,  13.0,    8.0},
        { 0.001,  0.001, 0.15,   0.1,    0.8,   8.0,   8.0,    5.0},
        { 0.001,  0.001, 0.15,   0.1,    0.8,   8.0,   8.0,    5.0},
        { 0.001,  0.001, 0.05,   0.22,   0.1,   0.16,  0.5,    0.9},
        { 0.001,  0.001, 0.05,   0.22,   0.1,   0.16,  0.5,    0.9}
    };

    if (ihole <= 0) {
        gamach = 0.0;
        std::ostringstream oss;
        oss << " No hole in SETGAM, gamach = " << gamach;
        logger().wlog(oss.str());
        return;
    }

    if (ihole > 16) {
        logger().wlog(" This version of FEFF will set gamach = 0.1 eV "
                      " for O1 and higher hole");
        logger().wlog(" You can use CORRECTIONS card to set "
                      " gamach = 0.1 + 2*vicorr");
    }

    double zz = static_cast<double>(iz);

    if (ihole <= 16) {
        // Build local arrays for terp interpolation (log10 of gamma vs Z)
        double zk[8], gamkp[8];
        int idx = ihole - 1;  // 0-based index into table
        for (int i = 0; i < 8; ++i) {
            gamkp[i] = std::log10(gamh[idx][i]);
            zk[i] = zh[idx][i];
        }
        feff::math::terp(zk, gamkp, 8, 1, zz, gamach);
    } else {
        // O-holes and higher: gamach = 0.1 eV (log10(0.1) = -1)
        gamach = -1.0;
    }

    // Convert from log10(gamma) to gamma
    gamach = std::pow(10.0, gamach);
}

void iniptz(FeffComplex ptz[3][3], int iptz, int modus) {
    // ptz is indexed [i+1][j+1] for i,j in {-1,0,1}
    // so ptz[-1] → ptz[0], ptz[0] → ptz[1], ptz[1] → ptz[2]

    constexpr FeffComplex zero(0.0, 0.0);
    constexpr FeffComplex one(1.0, 0.0);
    constexpr FeffComplex coni(0.0, 1.0);
    constexpr double sqrt2 = 1.4142135623730950488;

    if (iptz < 1 || iptz > 10) {
        logger().wlog("Inieln sees weird iptz - returns without "
                      "changing ptz - danger of calculating nonsense !!");
        return;
    }

    // Unity/3 matrix for orientation averaging
    double unity[3][3] = {};
    for (int i = 0; i < 3; ++i) {
        unity[i][i] = 1.0 / 3.0;
    }

    // Zero ptz
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            ptz[i][j] = zero;
        }
    }

    if (modus == 1) {
        // Spherical coordinates
        if (iptz == 10) {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    ptz[i][j] = FeffComplex(unity[i][j], 0.0);
                }
            }
        } else {
            // Map iptz (1-9) to row,col in -1..1 range
            int row = (iptz - 1) / 3;      // 0,1,2 → maps to -1,0,1
            int col = (iptz - 1) % 3;      // 0,1,2 → maps to -1,0,1
            ptz[row][col] = one;
        }
    } else if (modus == 2) {
        // Cartesian coordinates
        // Notation: ptz[i+1][j+1] where i,j in {-1,0,1}
        // idx: -1→0, 0→1, 1→2
        if (iptz == 10) {
            // Orientation average
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    ptz[i][j] = FeffComplex(unity[i][j], 0.0);
                }
            }
        } else if (iptz == 1) {
            // x x*
            ptz[2][2] = one / 2.0;    // (1,1)
            ptz[0][0] = one / 2.0;    // (-1,-1)
            ptz[0][2] = -one / 2.0;   // (-1,1)
            ptz[2][0] = -one / 2.0;   // (1,-1)
        } else if (iptz == 5) {
            // y y*
            ptz[2][2] = one / 2.0;    // (1,1)
            ptz[0][0] = one / 2.0;    // (-1,-1)
            ptz[0][2] = one / 2.0;    // (-1,1)
            ptz[2][0] = one / 2.0;    // (1,-1)
        } else if (iptz == 9) {
            // z z*
            ptz[1][1] = one;          // (0,0)
        } else if (iptz == 2) {
            // x y*
            ptz[2][2] = one * coni / 2.0;     // (1,1)
            ptz[0][0] = -one * coni / 2.0;    // (-1,-1)
            ptz[0][2] = -one * coni / 2.0;    // (-1,1)
            ptz[2][0] = one * coni / 2.0;     // (1,-1)
        } else if (iptz == 4) {
            // x* y
            ptz[2][2] = -one * coni / 2.0;    // (1,1)
            ptz[0][0] = one * coni / 2.0;     // (-1,-1)
            ptz[0][2] = -one * coni / 2.0;    // (-1,1)
            ptz[2][0] = one * coni / 2.0;     // (1,-1)
        } else if (iptz == 3) {
            // x z*
            ptz[0][1] = one / sqrt2;          // (-1,0)
            ptz[2][1] = -one / sqrt2;         // (1,0)
        } else if (iptz == 7) {
            // x* z
            ptz[1][0] = one / sqrt2;          // (0,-1)
            ptz[1][2] = -one / sqrt2;         // (0,1)
        } else if (iptz == 6) {
            // y z*
            ptz[0][1] = -one * coni / sqrt2;  // (-1,0)
            ptz[2][1] = -one * coni / sqrt2;  // (1,0)
        } else if (iptz == 8) {
            // y* z
            ptz[1][0] = one * coni / sqrt2;   // (0,-1)
            ptz[1][2] = one * coni / sqrt2;   // (0,1)
        }
    } else {
        throw std::runtime_error("alien modus in iniptz");
    }
}

} // namespace feff::common
