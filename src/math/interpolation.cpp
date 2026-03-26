#include "interpolation.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace feff::math {

int locat(double x, int n, const double xx[]) {
    int lo = 0;
    int u = n + 1;  // Must match Fortran: u = n+1 (not n)

    while (u - lo > 1) {
        int m = (u + lo) / 2;
        if (x < xx[m - 1]) {  // Fortran 1-based: xx(m) → xx[m-1]
            u = m;
        } else {
            lo = m;
        }
    }
    return lo;  // returns 0..n matching Fortran convention
}

void polint(const double xa[], const double ya[], int n, double x,
            double& y, double& dy) {
    constexpr int nmax = 4;
    double c[nmax], d[nmax];

    int ns = 0;
    double dif = std::abs(x - xa[0]);
    for (int i = 0; i < n; ++i) {
        double dift = std::abs(x - xa[i]);
        if (dift < dif) {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }
    y = ya[ns];
    ns--;

    for (int m = 1; m < n; ++m) {
        for (int i = 0; i < n - m; ++i) {
            double ho = xa[i] - x;
            double hp = xa[i + m] - x;
            double w = c[i + 1] - d[i];
            double den = ho - hp;
            if (den == 0.0) {
                throw std::runtime_error("failure in polint");
            }
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        if (2 * (ns + 1) < n - m) {
            dy = c[ns + 1];
        } else {
            dy = d[ns];
            ns--;
        }
        y += dy;
    }
}

void polinc(const double xa[], const FeffComplex ya[], int n, double x,
            FeffComplex& y, FeffComplex& dy) {
    constexpr int nmax = 4;
    FeffComplex c[nmax], d[nmax];

    int ns = 0;
    double dif = std::abs(x - xa[0]);
    for (int i = 0; i < n; ++i) {
        double dift = std::abs(x - xa[i]);
        if (dift < dif) {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }
    y = ya[ns];
    ns--;

    for (int m = 1; m < n; ++m) {
        for (int i = 0; i < n - m; ++i) {
            double ho = xa[i] - x;
            double hp = xa[i + m] - x;
            FeffComplex w = c[i + 1] - d[i];
            double den_d = ho - hp;
            if (den_d == 0.0) {
                throw std::runtime_error("failure in polinc");
            }
            FeffComplex den = w / den_d;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        if (2 * (ns + 1) < n - m) {
            dy = c[ns + 1];
        } else {
            dy = d[ns];
            ns--;
        }
        y += dy;
    }
}

void terp(const double x[], const double y[], int n, int m, double x0, double& y0) {
    // locat returns a 1-based index (matching Fortran convention):
    //   0 means x < x[0], i means x is between x[i-1] and x[i] (0-based),
    //   n means x >= x[n-1].
    // Fortran terp: k = min(max(i - m/2, 1), n - m) then polint(x(k),...) [1-based]
    // C++ equivalent: convert to 0-based by subtracting 1 from the Fortran k
    int i = locat(x0, n, x);
    int k = std::min(std::max(i - m / 2, 1), n - m) - 1;  // 0-based
    double dy;
    polint(x + k, y + k, m + 1, x0, y0, dy);
}

void terpc(const double x[], const FeffComplex y[], int n, int m,
           double x0, FeffComplex& y0) {
    int i = locat(x0, n, x);
    int k = std::min(std::max(i - m / 2, 1), n - m) - 1;  // 0-based
    FeffComplex dy;
    polinc(x + k, y + k, m + 1, x0, y0, dy);
}

int locat1(double x, int n, const float xx[]) {
    int lo = 0;
    int u = n + 1;  // Must match Fortran: u = n+1 (not n)

    while (u - lo > 1) {
        int m = (u + lo) / 2;
        if (x < xx[m - 1]) {
            u = m;
        } else {
            lo = m;
        }
    }
    return lo;
}

void terp1(const float x[], const float y[], int n, double x0, double& y0) {
    int i = locat1(x0, n, x);
    // Clamp to valid range
    i = std::max(i, 1);
    i = std::min(i, n - 1);

    // Convert to 0-based
    int idx = i - 1;

    if (x[idx + 1] - x[idx] == 0.0f) {
        throw std::runtime_error("TERP-1: division by zero");
    }

    y0 = y[idx] + (x0 - x[idx]) * (y[idx + 1] - y[idx]) / (x[idx + 1] - x[idx]);
}

} // namespace feff::math
