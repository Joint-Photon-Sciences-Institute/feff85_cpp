#include "sommerfeld.hpp"
#include <cmath>

namespace feff::math {

void somm(const double dr[], const double dp[], const double dq[],
          double dpas, double& da, int m, int np) {
    int mm = m + 1;
    double d1 = da + mm;
    da = 0.0;
    double db = 0.0;

    for (int i = 0; i < np; ++i) {
        double dl = std::pow(dr[i], mm);
        if (i != 0 && i != np - 1) {
            dl += dl;
            if (i % 2 == 1) dl += dl;  // Fortran: even i (1-based) → odd i (0-based)
        }
        double dc = dp[i] * dl;
        if (dc < 0.0) {
            db += dc;
        } else if (dc > 0.0) {
            da += dc;
        }
        dc = dq[i] * dl;
        if (dc < 0.0) {
            db += dc;
        } else if (dc > 0.0) {
            da += dc;
        }
    }
    da = dpas * (da + db) / 3.0;
    double dc_val = std::exp(dpas) - 1.0;
    db = d1 * (d1 + 1.0) * dc_val * std::exp((d1 - 1.0) * dpas);
    db = dr[0] * std::pow(dr[1], m) / db;
    dc_val = std::pow(dr[0], mm) * (1.0 + 1.0 / (dc_val * (d1 + 1.0))) / d1;
    da = da + dc_val * (dp[0] + dq[0]) - db * (dp[1] + dq[1]);
}

void somm2(const double dr[], const double dp[],
           double dpas, double& da, double rnrm, int m, int np) {
    int mm = m + 1;
    double d1 = da + mm;
    da = 0.0;

    double a1 = std::log(rnrm / dr[np - 3]) / dpas;
    double a2 = a1 * a1 / 8.0;
    double a3 = a1 * a1 * a1 / 12.0;

    for (int i = 0; i < np; ++i) {
        double dc;
        double drm = std::pow(dr[i], mm);
        if (i == 0) {
            dc = dp[i] * drm * 9.0 / 24.0;
        } else if (i == 1) {
            dc = dp[i] * drm * 28.0 / 24.0;
        } else if (i == 2) {
            dc = dp[i] * drm * 23.0 / 24.0;
        } else if (i == np - 4) {
            dc = dp[i] * drm * (25.0 / 24.0 - a2 + a3);
        } else if (i == np - 3) {
            dc = dp[i] * drm * (0.5 + a1 - 3.0 * a2 - a3);
        } else if (i == np - 2) {
            dc = dp[i] * drm * (-1.0 / 24.0 + 5.0 * a2 - a3);
        } else if (i == np - 1) {
            dc = dp[i] * drm * (-a2 + a3);
        } else {
            dc = dp[i] * drm;
        }
        da += dc;
    }
    da = dpas * da;

    // Initial point correction
    double dd = std::exp(dpas) - 1.0;
    double db = d1 * (d1 + 1.0) * dd * std::exp((d1 - 1.0) * dpas);
    db = dr[0] * std::pow(dr[1], m) / db;
    dd = std::pow(dr[0], mm) * (1.0 + 1.0 / (dd * (d1 + 1.0))) / d1;
    da = da + dd * dp[0] - db * dp[1];
}

void csomm(const double dr[], const FeffComplex dp[], const FeffComplex dq[],
           double dpas, FeffComplex& da, int m, int np) {
    int mm = m + 1;
    double d1 = da.real() + mm;
    da = FeffComplex(0.0, 0.0);

    for (int i = 0; i < np; ++i) {
        double dl = std::pow(dr[i], mm);
        if (i != 0 && i != np - 1) {
            dl += dl;
            if (i % 2 == 1) dl += dl;
        }
        FeffComplex dc = dp[i] * dl;
        da += dc;
        dc = dq[i] * dl;
        da += dc;
    }
    da = dpas * da / 3.0;
    double dd = std::exp(dpas) - 1.0;
    double db = d1 * (d1 + 1.0) * dd * std::exp((d1 - 1.0) * dpas);
    db = dr[0] * std::pow(dr[1], m) / db;
    dd = std::pow(dr[0], mm) * (1.0 + 1.0 / (dd * (d1 + 1.0))) / d1;
    da = da + dd * (dp[0] + dq[0]) - db * (dp[1] + dq[1]);
}

void csomm2(const double dr[], const FeffComplex dp[],
            double dpas, FeffComplex& da, double rnrm, int np) {
    double d1 = da.real() + 1.0;
    da = FeffComplex(0.0, 0.0);

    double a1 = std::log(rnrm / dr[np - 3]) / dpas;
    double a2 = a1 * a1 / 8.0;
    double a3 = a1 * a1 * a1 / 12.0;

    for (int i = 0; i < np; ++i) {
        FeffComplex dc;
        if (i == 0) {
            dc = dp[i] * dr[i] * 9.0 / 24.0;
        } else if (i == 1) {
            dc = dp[i] * dr[i] * 28.0 / 24.0;
        } else if (i == 2) {
            dc = dp[i] * dr[i] * 23.0 / 24.0;
        } else if (i == np - 4) {
            dc = dp[i] * dr[i] * (25.0 / 24.0 - a2 + a3);
        } else if (i == np - 3) {
            dc = dp[i] * dr[i] * (0.5 + a1 - 3.0 * a2 - a3);
        } else if (i == np - 2) {
            dc = dp[i] * dr[i] * (-1.0 / 24.0 + 5.0 * a2 - a3);
        } else if (i == np - 1) {
            dc = dp[i] * dr[i] * (-a2 + a3);
        } else {
            dc = dp[i] * dr[i];
        }
        da += dc;
    }
    da = dpas * da;

    // Initial point correction
    double dd = std::exp(dpas) - 1.0;
    double db = d1 * (d1 + 1.0) * dd * std::exp((d1 - 1.0) * dpas);
    db = dr[0] / db;
    dd = dr[0] * (1.0 + 1.0 / (dd * (d1 + 1.0))) / d1;
    da = da + dd * dp[0] - db * dp[1];
}

} // namespace feff::math
