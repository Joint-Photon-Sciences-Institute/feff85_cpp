// Low-level numerical utilities for the ATOM module.
// Converted from: aprdev.f, dentfa.f, cofcon.f, potslw.f, messer.f, bkmrdf.f

#include "utility.hpp"
#include "../common/logging.hpp"
#include "../math/wigner.hpp"
#include <cmath>
#include <sstream>
#include <algorithm>

namespace feff::atom {

// =========================================================================
// aprdev — polynomial product coefficient
// Fortran: aprdev = sum_{m=1}^{l} a(m)*b(l+1-m)
// C++: arrays are 0-based, l is the Fortran 1-based index.
// =========================================================================
double aprdev(const double a[], const double b[], int l) {
    double result = 0.0;
    // Fortran: do m=1,l => aprdev = aprdev + a(m)*b(l+1-m)
    // C++ 0-based: m goes 0..l-1, indices a[m], b[l-1-m]
    for (int m = 0; m < l; ++m) {
        result += a[m] * b[l - 1 - m];
    }
    return result;
}

// =========================================================================
// dentfa — Thomas-Fermi potential approximation
// =========================================================================
double dentfa(double dr, double dz, double ch) {
    if ((dz + ch) < 1.0e-04) return 0.0;

    // Match Fortran single-precision constants exactly.
    // In Fortran dentfa.f, ALL literal constants (0.60112, 1.81061, 0.8853, etc.)
    // are single-precision (no 'd0' suffix). When used in mixed-mode expressions
    // with double-precision variables, they are first promoted from float to double,
    // which preserves only ~7 significant digits. Using C++ double literals would
    // give ~16-digit precision, producing a ~3e-8 relative difference that
    // propagates through the SCF loop and causes convergence to a wrong Fermi energy.
    double w = dr * std::pow(dz + ch, static_cast<double>(1.0f / 3.0f));
    w = std::sqrt(w / static_cast<double>(0.8853f));
    double t = w * (static_cast<double>(0.60112f) * w + static_cast<double>(1.81061f)) + 1.0;
    w = w * (w * (w * (w * (static_cast<double>(0.04793f) * w
        + static_cast<double>(0.21465f)) + static_cast<double>(0.77112f))
        + static_cast<double>(1.39515f)) + static_cast<double>(1.81061f)) + 1.0;
    return (dz + ch) * (1.0 - (t / w) * (t / w)) / dr;
}

// =========================================================================
// cofcon — convergence acceleration
// =========================================================================
void cofcon(double& a, double& b, double& p, double& q) {
    if (p * q < 0.0) {
        if (b >= 0.2) b -= 0.1;
    } else if (p * q > 0.0) {
        if (b <= 0.8) b += 0.1;
    }
    a = 1.0 - b;
    q = p;
}

// =========================================================================
// potslw — Coulomb potential from density (4-point integration)
// Arrays are 0-based in C++; Fortran indices [1..np] -> C++ [0..np-1].
// =========================================================================
void potslw(double dv[], const double d[], const double dr[], double dpas, int np) {
    // dp is a local work array (Fortran dimension 251)
    double dp[atom_grid];

    double das = dpas / 24.0;

    // dv(i) = d(i)*dr(i)  (temporary storage of r*rho)
    for (int i = 0; i < np; ++i) {
        dv[i] = d[i] * dr[i];
    }

    double dlo = std::exp(dpas);
    double dlo2 = dlo * dlo;

    // Fortran: dp(2) = dr(1)*(d(2)-d(1)*dlo2)/(12*(dlo-1))
    // C++ 0-based: dp[1] = dr[0]*(d[1]-d[0]*dlo2)/(12*(dlo-1))
    dp[1] = dr[0] * (d[1] - d[0] * dlo2) / (12.0 * (dlo - 1.0));
    dp[0] = dv[0] / 3.0 - dp[1] / dlo2;
    dp[1] = dv[1] / 3.0 - dp[1] * dlo2;

    // j = np-1 in Fortran (1-based) => j = np-2 in C++ (0-based)
    int j = np - 2;
    // Fortran: do i=3,j => C++: i from 2 to j-1
    for (int i = 2; i <= j; ++i) {
        dp[i] = dp[i - 1] + das * (13.0 * (dv[i] + dv[i - 1]) - (dv[i - 2] + dv[i + 1]));
    }

    dp[np - 1] = dp[j];
    dv[j] = dp[j];
    dv[np - 1] = dp[j];

    // Fortran: do i=3,j => k=np+1-i goes from np-2 down to 2 (1-based)
    // C++ 0-based: k goes from np-3 down to 1
    for (int k = np - 3; k >= 1; --k) {
        dv[k] = dv[k + 1] / dlo + das *
            (13.0 * (dp[k + 1] / dlo + dp[k]) - (dp[k + 2] / dlo2 + dp[k - 1] * dlo));
    }

    // Fortran: dv(1)=dv(3)/dlo2+dpas*(dp(1)+4*dp(2)/dlo+dp(3)/dlo2)/3
    dv[0] = dv[2] / dlo2 + dpas * (dp[0] + 4.0 * dp[1] / dlo + dp[2] / dlo2) / 3.0;

    for (int i = 0; i < np; ++i) {
        dv[i] = dv[i] / dr[i];
    }
}

// =========================================================================
// messer — error message printer
// =========================================================================
void messer(ErrorState& err) {
    int ilig = err.numerr / 1000;
    int ier = err.numerr - 1000 * ilig;

    std::ostringstream oss;
    oss << "error number " << ier
        << " detected on a line " << ilig
        << "in the program" << std::string(err.dlabpr, 8);
    feff::common::logger().wlog(oss.str());
}

// =========================================================================
// bkmrdf — Breit interaction angular coefficients
// Fortran i,j are 1-based orbital indices; C++ uses 0-based.
// =========================================================================
void bkmrdf(int i, int j, int k, OrbitalConfig& config, BreitCoefficients& breit) {
    for (int l = 0; l < 3; ++l) {
        breit.cmag[l] = 0.0;
        breit.cret[l] = 0.0;
    }

    int ji = 2 * std::abs(config.kap[i]) - 1;
    int jj = 2 * std::abs(config.kap[j]) - 1;
    int kam = config.kap[j] - config.kap[i];

    int l = k - 1;
    for (int m = 0; m < 3; ++m) {
        if (l >= 0) {
            double a = feff::math::cwig3j(ji, jj, l + l, -1, 1, 2);
            a = a * a;
            if (a != 0.0) {
                double c = static_cast<double>(l + l + 1);
                double cm, cz, cp, d;
                if (m == 1) {
                    // Fortran m=2 (middle case)
                    d = static_cast<double>(k * (k + 1));
                    cm = static_cast<double>((config.kap[i] + config.kap[j]) * (config.kap[i] + config.kap[j]));
                    cz = cm;
                    cp = cm;
                } else {
                    int n;
                    if (m == 0) {
                        // Fortran m=1 (m < 2)
                        cm = static_cast<double>((kam + k) * (kam + k));
                        cz = static_cast<double>(kam * kam - k * k);
                        cp = static_cast<double>((k - kam) * (k - kam));
                        n = k;
                    } else {
                        // m == 2 => Fortran m=3 (m > 2)
                        cm = static_cast<double>((kam - l) * (kam - l));
                        cz = static_cast<double>(kam * kam - l * l);
                        cp = static_cast<double>((kam + l) * (kam + l));
                        n = l;
                        c = -c;
                    }
                    int l1 = l + 1;
                    double am = static_cast<double>((kam - l) * (kam + l1)) / c;
                    double az = static_cast<double>(kam * kam + l * l1) / c;
                    double ap = static_cast<double>((l + kam) * (kam - l1)) / c;
                    d = static_cast<double>(n * (k + k + 1));
                    c = std::abs(c) * d;
                    if (c != 0.0) c = static_cast<double>(n) / c;
                    breit.cret[0] += a * (am - c * cm);
                    breit.cret[1] += (a + a) * (az - c * cz);
                    breit.cret[2] += a * (ap - c * cp);

                    // d was set above, fall through to cmag update
                }
                if (d != 0.0) {
                    double ad = a / d;
                    breit.cmag[0] += cm * ad;
                    breit.cmag[1] += cz * (ad + ad);
                    breit.cmag[2] += cp * ad;
                }
            }
        }
        l = l + 1;
    }
}

} // namespace feff::atom
