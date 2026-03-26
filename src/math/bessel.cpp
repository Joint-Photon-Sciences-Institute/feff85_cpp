#include "bessel.hpp"
#include <cmath>
#include <stdexcept>

namespace feff::math {

void bjnser(FeffComplex x, int l, FeffComplex& jl, FeffComplex& nl, int ifl) {
    constexpr int niter = 160;
    constexpr double tol = 1.0e-15;

    if (l < 0) {
        throw std::runtime_error("l < 0 in bjnser");
    }
    if (x.real() < 0.0) {
        throw std::runtime_error("Re(x) < 0 in bjnser");
    }

    int lp1 = l + 1;
    FeffComplex u = x * x / 2.0;

    // djl = 1 * 3 * 5 * ... * (2*l+1)
    double djl = 1.0;
    double fac = -1.0;
    for (int il = 0; il < lp1; ++il) {
        fac += 2.0;
        djl = fac * djl;
    }
    double dnl = djl / (2 * l + 1);

    // Calculate jl
    if (ifl != 2) {
        FeffComplex pj(1.0, 0.0);
        int nf = 1;
        int nfac = 2 * l + 3;
        double den = nfac;
        double sgn = -1.0;
        FeffComplex ux = u;
        bool converged = false;
        for (int il = 0; il < niter; ++il) {
            FeffComplex del = sgn * ux / den;
            pj += del;
            double trel = std::abs(del / pj);
            if (trel <= tol) {
                converged = true;
                break;
            }
            sgn = -sgn;
            ux *= u;
            nf += 1;
            nfac += 2;
            den = nf * nfac * den;
        }
        if (!converged) {
            throw std::runtime_error("jl does not converge in bjnser");
        }
        jl = pj * std::pow(x, l) / djl;
    }

    if (ifl == 1) return;

    // Calculate nl
    {
        FeffComplex pn(1.0, 0.0);
        int nf = 1;
        int nfac = 1 - 2 * l;
        double den = nfac;
        double sgn = -1.0;
        FeffComplex ux = u;
        bool converged = false;
        for (int il = 0; il < niter; ++il) {
            FeffComplex del = sgn * ux / den;
            pn += del;
            double trel = std::abs(del / pn);
            if (trel <= tol) {
                converged = true;
                break;
            }
            sgn = -sgn;
            ux *= u;
            nf += 1;
            nfac += 2;
            den = nf * nfac * den;
        }
        if (!converged) {
            throw std::runtime_error("nl does not converge in bjnser");
        }
        nl = -pn * dnl / std::pow(x, lp1);
    }
}

void besjn(FeffComplex x, FeffComplex jl[], FeffComplex nl[]) {
    constexpr double xcut = 1.0;
    constexpr double xcut1 = 7.51;
    constexpr double xcut2 = 5.01;

    if (x.real() <= 0.0) {
        throw std::runtime_error("Re(x) <= 0 in besjn");
    }

    int lmaxp1 = ltot + 2;

    if (x.real() < xcut && std::abs(x.imag()) < xcut) {
        // Case Re(x) < 1: use series expansion
        for (int il = 0; il < lmaxp1; ++il) {
            FeffComplex xjl, xnl;
            bjnser(x, il, xjl, xnl, 0);
            jl[il] = xjl;
            nl[il] = xnl;
        }
    } else if (x.real() < xcut1 && std::abs(x.imag()) < xcut1) {
        // Case 1 <= Re(x) < 7.5
        FeffComplex xjl, xnl;
        bjnser(x, lmaxp1 - 1, xjl, xnl, 1);
        jl[lmaxp1 - 1] = xjl;

        bjnser(x, lmaxp1 - 2, xjl, xnl, 1);
        jl[lmaxp1 - 2] = xjl;

        if (x.real() < xcut2 && std::abs(x.imag()) < xcut2) {
            bjnser(x, 0, xjl, xnl, 2);
            nl[0] = xnl;
            bjnser(x, 1, xjl, xnl, 2);
            nl[1] = xnl;
        } else {
            FeffComplex asx = std::sin(x);
            FeffComplex acx = std::cos(x);
            FeffComplex xi = 1.0 / x;
            FeffComplex xi2 = xi * xi;
            nl[0] = -acx * xi;
            nl[1] = -acx * xi2 - asx * xi;
        }

        // Recursion for nl (forward) - Fortran: nl(lp1) for lp1=3..lmaxp1
        // In 0-based: nl[i] for i=2..lmaxp1-1
        for (int i = 2; i < lmaxp1; ++i) {
            int l = i - 1;  // Fortran l = lp1 - 2
            double tlxp1 = 2 * l + 1;
            nl[i] = tlxp1 * nl[i - 1] / x - nl[i - 2];
        }

        // Recursion for jl (backward) - Fortran: jl(lp1) for lp1=lmaxp1-2 down to 1
        // In 0-based: jl[i] for i=lmaxp1-3 down to 0
        for (int i = lmaxp1 - 3; i >= 0; --i) {
            int l = i;  // l value for this element
            double tlxp3 = 2 * (l + 1) + 1;  // 2*l+3 where l=lp1-1 in Fortran
            jl[i] = tlxp3 * jl[i + 1] / x - jl[i + 2];
        }
    } else {
        // Case Re(x) > 7.5: asymptotic expansion
        FeffComplex sjl[ltot + 2], cjl_arr[ltot + 2], snl_arr[ltot + 2], cnl_arr[ltot + 2];

        FeffComplex xi = 1.0 / x;
        FeffComplex xi2 = xi * xi;
        FeffComplex xi3 = xi * xi2;
        FeffComplex xi4 = xi * xi3;
        FeffComplex xi5 = xi * xi4;
        FeffComplex xi6 = xi * xi5;
        FeffComplex xi7 = xi * xi6;
        FeffComplex xi8 = xi * xi7;
        FeffComplex xi9 = xi * xi8;
        FeffComplex xi10 = xi * xi9;
        FeffComplex xi11 = xi * xi10;

        sjl[0] = xi;
        sjl[1] = xi2;
        sjl[2] = 3.0 * xi3 - xi;
        sjl[3] = 15.0 * xi4 - 6.0 * xi2;
        sjl[4] = 105.0 * xi5 - 45.0 * xi3 + xi;
        sjl[5] = 945.0 * xi6 - 420.0 * xi4 + 15.0 * xi2;
        sjl[6] = 10395.0 * xi7 - 4725.0 * xi5 + 210.0 * xi3 - xi;
        sjl[7] = 135135.0 * xi8 - 62370.0 * xi6 + 3150.0 * xi4 - 28.0 * xi2;
        sjl[8] = 2027025.0 * xi9 - 945945.0 * xi7 + 51975.0 * xi5
                 - 630.0 * xi3 + xi;
        sjl[9] = 34459425.0 * xi10 - 16216200.0 * xi8 + 945945.0 * xi6
                 - 13860.0 * xi4 + 45.0 * xi2;
        sjl[10] = 654729075.0 * xi11 - 310134825.0 * xi9 + 18918900.0 * xi7
                  - 315315.0 * xi5 + 1485.0 * xi3 - xi;

        cjl_arr[0] = FeffComplex(0.0, 0.0);
        cjl_arr[1] = -xi;
        cjl_arr[2] = -3.0 * xi2;
        cjl_arr[3] = -15.0 * xi3 + xi;
        cjl_arr[4] = -105.0 * xi4 + 10.0 * xi2;
        cjl_arr[5] = -945.0 * xi5 + 105.0 * xi3 - xi;
        cjl_arr[6] = -10395.0 * xi6 + 1260.0 * xi4 - 21.0 * xi2;
        cjl_arr[7] = -135135.0 * xi7 + 17325.0 * xi5 - 378.0 * xi3 + xi;
        cjl_arr[8] = -2027025.0 * xi8 + 270270.0 * xi6 - 6930.0 * xi4 + 36.0 * xi2;
        cjl_arr[9] = -34459425.0 * xi9 + 4729725.0 * xi7 - 135135.0 * xi5
                     + 990.0 * xi3 - xi;
        cjl_arr[10] = -654729075.0 * xi10 + 91891800.0 * xi8 - 2837835.0 * xi6
                      + 25740.0 * xi4 - 55.0 * xi2;

        for (int ie = 0; ie < 11; ++ie) {
            snl_arr[ie] = cjl_arr[ie];
            cnl_arr[ie] = -sjl[ie];
        }

        // Recursion for higher orders
        for (int lp1 = 11; lp1 < lmaxp1; ++lp1) {
            int l = lp1 - 1;
            double tlxp1 = static_cast<double>(2 * l + 1);
            sjl[lp1] = tlxp1 * xi * sjl[lp1 - 1] - sjl[lp1 - 2];
            cjl_arr[lp1] = tlxp1 * xi * cjl_arr[lp1 - 1] - cjl_arr[lp1 - 2];
            snl_arr[lp1] = tlxp1 * xi * snl_arr[lp1 - 1] - snl_arr[lp1 - 2];
            cnl_arr[lp1] = tlxp1 * xi * cnl_arr[lp1 - 1] - cnl_arr[lp1 - 2];
        }

        FeffComplex asx = std::sin(x);
        FeffComplex acx = std::cos(x);
        for (int lp1 = 0; lp1 < lmaxp1; ++lp1) {
            jl[lp1] = asx * sjl[lp1] + acx * cjl_arr[lp1];
            nl[lp1] = asx * snl_arr[lp1] + acx * cnl_arr[lp1];
        }
    }
}

void besjh(FeffComplex x, int lbmax, FeffComplex jl[], FeffComplex hl[]) {
    constexpr double xcut = 1.0;
    constexpr double xcut1 = 7.51;
    constexpr double xcut2 = 5.01;

    if (x.real() < 0.0) {
        throw std::runtime_error("Re(x) < 0 in besjh");
    }

    int lmax = std::min(lbmax, ltot + 1);
    int lmaxp1 = lmax + 1;

    // Temporary nl array for intermediate computation
    FeffComplex nl_tmp[ltot + 2];

    if (x.real() < xcut && std::abs(x.imag()) < xcut) {
        // Case |Re(x)| < 1 and |Im(x)| < 1: series expansion for all l
        for (int ll = 0; ll <= lmax; ++ll) {
            FeffComplex xjl, xnl;
            bjnser(x, ll, xjl, xnl, 0);
            jl[ll] = xjl;
            hl[ll] = -xnl + coni * xjl;
        }
    } else if (x.real() < xcut && std::abs(x.imag()) < xcut1) {
        // Case Re(x) < 1 but |Im(x)| >= 1: use direct series for all l
        // to avoid unstable backward recursion for jl with large imaginary args
        for (int ll = 0; ll <= lmax; ++ll) {
            FeffComplex xjl, xnl;
            bjnser(x, ll, xjl, xnl, 0);
            jl[ll] = xjl;
            hl[ll] = -xnl + coni * xjl;
        }
    } else if (x.real() < xcut1 && std::abs(x.imag()) < xcut1) {
        // Case 1 <= Re(x) < 7.5
        FeffComplex xjl, xnl;
        bjnser(x, lmax, xjl, xnl, 1);
        jl[lmax] = xjl;

        bjnser(x, lmax - 1, xjl, xnl, 1);
        jl[lmax - 1] = xjl;

        if (x.real() < xcut2 && std::abs(x.imag()) < xcut2) {
            bjnser(x, 0, xjl, xnl, 2);
            nl_tmp[0] = xnl;
            bjnser(x, 1, xjl, xnl, 2);
            nl_tmp[1] = xnl;
        } else {
            FeffComplex asx = std::sin(x);
            FeffComplex acx = std::cos(x);
            FeffComplex xi = 1.0 / x;
            FeffComplex xi2 = xi * xi;
            nl_tmp[0] = -acx * xi;
            nl_tmp[1] = -acx * xi2 - asx * xi;
        }

        // Forward recursion for nl
        for (int lp1 = 2; lp1 < lmaxp1; ++lp1) {
            int l = lp1 - 1;
            double tlxp1 = 2 * l + 1;
            nl_tmp[lp1] = tlxp1 * nl_tmp[lp1 - 1] / x - nl_tmp[lp1 - 2];
        }

        // Backward recursion for jl
        // Fortran: do lxx=3,lmaxp1; lp1=lmaxp1+1-lxx; l=lp1-1; jl(l)=(2l+3)*jl(l+1)/x-jl(l+2)
        for (int lxx = 2; lxx < lmaxp1; ++lxx) {
            int lp1 = lmaxp1 - lxx;
            int l = lp1 - 1;  // FIX: was lp1 (off by one)
            double tlxp3 = 2 * l + 3;
            jl[l] = tlxp3 * jl[l + 1] / x - jl[l + 2];
        }

        // Compute hl from jl and nl
        for (int il = 0; il < lmaxp1; ++il) {
            hl[il] = -nl_tmp[il] + coni * jl[il];
        }
    } else {
        // Case Re(x) > 7.5: asymptotic expansion
        FeffComplex sjl[ltot + 2], cjl_arr[ltot + 2];

        FeffComplex xi = 1.0 / x;
        FeffComplex xi2 = xi * xi;
        FeffComplex xi3 = xi * xi2;
        FeffComplex xi4 = xi * xi3;
        FeffComplex xi5 = xi * xi4;
        FeffComplex xi6 = xi * xi5;
        FeffComplex xi7 = xi * xi6;
        FeffComplex xi8 = xi * xi7;
        FeffComplex xi9 = xi * xi8;
        FeffComplex xi10 = xi * xi9;
        FeffComplex xi11 = xi * xi10;

        sjl[0] = xi;
        sjl[1] = xi2;
        sjl[2] = 3.0 * xi3 - xi;
        sjl[3] = 15.0 * xi4 - 6.0 * xi2;
        sjl[4] = 105.0 * xi5 - 45.0 * xi3 + xi;
        sjl[5] = 945.0 * xi6 - 420.0 * xi4 + 15.0 * xi2;
        sjl[6] = 10395.0 * xi7 - 4725.0 * xi5 + 210.0 * xi3 - xi;
        sjl[7] = 135135.0 * xi8 - 62370.0 * xi6 + 3150.0 * xi4 - 28.0 * xi2;
        sjl[8] = 2027025.0 * xi9 - 945945.0 * xi7 + 51975.0 * xi5
                 - 630.0 * xi3 + xi;
        sjl[9] = 34459425.0 * xi10 - 16216200.0 * xi8 + 945945.0 * xi6
                 - 13860.0 * xi4 + 45.0 * xi2;
        sjl[10] = 654729075.0 * xi11 - 310134825.0 * xi9 + 18918900.0 * xi7
                  - 315315.0 * xi5 + 1485.0 * xi3 - xi;

        cjl_arr[0] = FeffComplex(0.0, 0.0);
        cjl_arr[1] = -xi;
        cjl_arr[2] = -3.0 * xi2;
        cjl_arr[3] = -15.0 * xi3 + xi;
        cjl_arr[4] = -105.0 * xi4 + 10.0 * xi2;
        cjl_arr[5] = -945.0 * xi5 + 105.0 * xi3 - xi;
        cjl_arr[6] = -10395.0 * xi6 + 1260.0 * xi4 - 21.0 * xi2;
        cjl_arr[7] = -135135.0 * xi7 + 17325.0 * xi5 - 378.0 * xi3 + xi;
        cjl_arr[8] = -2027025.0 * xi8 + 270270.0 * xi6 - 6930.0 * xi4 + 36.0 * xi2;
        cjl_arr[9] = -34459425.0 * xi9 + 4729725.0 * xi7 - 135135.0 * xi5
                     + 990.0 * xi3 - xi;
        cjl_arr[10] = -654729075.0 * xi10 + 91891800.0 * xi8 - 2837835.0 * xi6
                      + 25740.0 * xi4 - 55.0 * xi2;

        // Recursion for higher orders
        for (int lp1 = 11; lp1 < lmaxp1; ++lp1) {
            int l = lp1 - 1;
            double tlxp1 = static_cast<double>(2 * l + 1);
            sjl[lp1] = tlxp1 * xi * sjl[lp1 - 1] - sjl[lp1 - 2];
            cjl_arr[lp1] = tlxp1 * xi * cjl_arr[lp1 - 1] - cjl_arr[lp1 - 2];
        }

        FeffComplex asx = std::sin(x);
        FeffComplex acx = std::cos(x);
        FeffComplex epx;
        if (x.imag() >= 0.0) {
            epx = std::exp(coni * x);
        } else {
            epx = std::exp(-coni * x);
        }

        for (int ll = 0; ll <= lmax; ++ll) {
            jl[ll] = asx * sjl[ll] + acx * cjl_arr[ll];
            if (x.imag() >= 0.0) {
                hl[ll] = (sjl[ll] + coni * cjl_arr[ll]) * epx;
            } else {
                hl[ll] = (sjl[ll] - coni * cjl_arr[ll]) * epx;
            }
        }
    }
}

void exjlnl(FeffComplex z, int l, FeffComplex& jl, FeffComplex& nl) {
    // For small |z|, use series expansion
    if (std::abs(z) < 0.3) {
        bjnser(z, l, jl, nl, 0);
        return;
    }

    FeffComplex cosz = std::cos(z);
    FeffComplex sinz = std::sin(z);

    if (l == 0) {
        jl = sinz / z;
        nl = -cosz / z;
    } else if (l == 1) {
        jl = sinz / (z * z) - cosz / z;
        nl = -cosz / (z * z) - sinz / z;
    } else if (l == 2) {
        FeffComplex z2 = z * z;
        FeffComplex z3 = z2 * z;
        jl = (3.0 / z3 - 1.0 / z) * sinz - 3.0 * cosz / z2;
        nl = (-3.0 / z3 + 1.0 / z) * cosz - 3.0 * sinz / z2;
    } else if (l == 3) {
        FeffComplex z2 = z * z;
        FeffComplex z3 = z2 * z;
        FeffComplex z4 = z3 * z;
        jl = (15.0 / z4 - 6.0 / z2) * sinz + (-15.0 / z3 + 1.0 / z) * cosz;
        nl = (-15.0 / z4 + 6.0 / z2) * cosz + (-15.0 / z3 + 1.0 / z) * sinz;
    } else if (l == 4) {
        FeffComplex z2 = z * z;
        FeffComplex z3 = z2 * z;
        FeffComplex z4 = z3 * z;
        FeffComplex z5 = z4 * z;
        jl = (105.0 / z5 - 45.0 / z3 + 1.0 / z) * sinz +
             (-105.0 / z4 + 10.0 / z2) * cosz;
        nl = (-105.0 / z5 + 45.0 / z3 - 1.0 / z) * cosz +
             (-105.0 / z4 + 10.0 / z2) * sinz;
    } else if (l == 5) {
        FeffComplex z2 = z * z;
        FeffComplex z3 = z2 * z;
        FeffComplex z4 = z3 * z;
        FeffComplex z5 = z4 * z;
        FeffComplex z6 = z5 * z;
        jl = (945.0 / z6 - 420.0 / z4 + 15.0 / z2) * sinz +
             (-945.0 / z5 + 105.0 / z3 - 1.0 / z) * cosz;
        nl = (-945.0 / z6 + 420.0 / z4 - 15.0 / z2) * cosz +
             (-945.0 / z5 + 105.0 / z3 - 1.0 / z) * sinz;
    } else if (l == 6) {
        FeffComplex z2 = z * z;
        FeffComplex z3 = z2 * z;
        FeffComplex z4 = z3 * z;
        FeffComplex z5 = z4 * z;
        FeffComplex z6 = z5 * z;
        FeffComplex z7 = z6 * z;
        jl = (10395.0 / z7 - 4725.0 / z5 + 210.0 / z3 - 1.0 / z) * sinz +
             (-10395.0 / z6 + 1260.0 / z4 - 21.0 / z2) * cosz;
        nl = (-10395.0 / z7 + 4725.0 / z5 - 210.0 / z3 + 1.0 / z) * cosz +
             (-10395.0 / z6 + 1260.0 / z4 - 21.0 / z2) * sinz;
    } else {
        throw std::runtime_error("exjlnl: l out of range (must be 0-6)");
    }
}

} // namespace feff::math
