// Nuclear potential construction and radial mesh setup.
// Converted from: nucdev.f

#include "nucdev.hpp"
#include "../par/parallel.hpp"
#include <cmath>

namespace feff::atom {

// Integer power matching Fortran x**n behavior
static inline double fpow_int(double x, int n) {
    if (n == 0) return 1.0;
    if (n < 0) { x = 1.0/x; n = -n; }
    double r = 1.0;
    for (int i = 0; i < n; ++i) r *= x;
    return r;
}

void nucdev(double av[], double dr[], double dv[], double dz, double hx,
            int& nuc, int np, int ndor, double& dr1) {
    // Local work array (Fortran: at(251))
    double at[atom_grid];

    // Specify atomic mass and thickness of nuclear shell.
    // a = atomic mass (<=0 for point charge)
    // epai = Fermi density parameter (<=0 for uniform distribution)
    double a = 0.0;
    double epai = 0.0;

    // Calculate nuclear radius index
    if (a <= 0.1) {
        nuc = 1;
    } else {
        // Finite nucleus case (currently dead code since a=0)
        // Fortran: a**(1./3.) — single-precision exponent (1./3. = 0.333333343267...)
        // Currently dead code (a=0 above), but match Fortran for correctness.
        a = dz * std::pow(a, static_cast<double>(1.0f / 3.0f)) * 2.2677e-05;
        // Fortran nuc is 1-based; here we keep the same convention internally
        // and convert at boundaries. nuc passed in is 1-based from caller.
        double b = a / std::exp(hx * (nuc - 1));
        if (b <= dr1) {
            dr1 = b;
        } else {
            b = std::log(a / dr1) / hx;
            nuc = 3 + 2 * static_cast<int>(b / 2.0);
            if (nuc >= np) {
                feff::par::par_stop("dr1 too small");
            }
            dr1 = a * std::exp(-(nuc - 1) * hx);
        }
    }

    // Build radial mesh: dr[i] = dr1/dz * exp(hx * i)  (0-based)
    // Fortran: dr(1) = dr1/dz, dr(l) = dr(1)*exp(hx*(l-1))
    dr[0] = dr1 / dz;
    for (int l = 1; l < np; ++l) {
        dr[l] = dr[0] * std::exp(hx * l);
    }

    if (ndor < 5) {
        feff::par::par_stop("stopped in programm nucdev, ndor should be > 4.");
    }

    // Zero out development coefficients
    for (int i = 0; i < ndor; ++i) {
        av[i] = 0.0;
    }

    if (epai <= 0.0) {
        // Uniform or point charge distribution
        for (int i = 0; i < np; ++i) {
            dv[i] = -dz / dr[i];
        }

        if (nuc <= 1) {
            // Point charge: av[0] = -dz (Fortran av(1)=-dz)
            av[0] = -dz;
        } else {
            // Uniform finite nucleus
            // Fortran: nuc is 1-based index into dr
            // av(2) = -3*dz/(2*dr(nuc))
            int nuc0 = nuc - 1;  // 0-based index
            av[1] = -3.0 * dz / (dr[nuc0] + dr[nuc0]);
            // av(4) = -av(2)/(3*dr(nuc)^2)
            av[3] = -av[1] / (3.0 * dr[nuc0] * dr[nuc0]);
            // Fortran: do i=1,nuc-1 => overwrite dv for interior points
            for (int i = 0; i < nuc0; ++i) {
                dv[i] = av[1] + av[3] * dr[i] * dr[i];
            }
        }
    } else {
        // Fermi distribution nuclear charge (epai > 0 case)
        // nuc0 = 0-based index of nuclear radius
        int nuc0 = nuc - 1;

        double b = std::exp(-dr[nuc0] / epai);
        b = 1.0 / (1.0 + b);
        av[3] = b;   // Fortran av(4)
        av[4] = epai * b * (b - 1.0);  // Fortran av(5)

        if (ndor > 5) {
            at[0] = 1.0;  // Fortran at(1)
            at[1] = 1.0;  // Fortran at(2)
            int nf = 1;
            // Fortran loop: do i=6,ndor (1-based) => C++ i=5..ndor-1 (0-based)
            for (int ii = 5; ii < ndor; ++ii) {
                // Fortran: n = i - 4  (where i is 1-based, so i_f = ii+1)
                int n = (ii + 1) - 4;  // = ii - 3
                nf = n * nf;
                int n1 = n + 1;

                // Use dv[] as scratch; all 1-based Fortran indexing
                // dv(1) = n * at(1)
                dv[0] = n * at[0];
                // dv(n1) = 1.0
                dv[n1 - 1] = 1.0;
                // do j=2,n: dv(j) = (n-j+2)*at(j-1) + (n-j+1)*at(j)
                for (int jf = 2; jf <= n; ++jf) {
                    dv[jf - 1] = (n - jf + 2) * at[jf - 2] + (n - jf + 1) * at[jf - 1];
                }

                // do j=1,n1: accumulate av(i) and copy dv->at
                av[ii] = 0.0;
                for (int jf = 1; jf <= n1; ++jf) {
                    int mm = n + 1 - jf;
                    int sgn = (jf % 2 == 0) ? -1 : 1;
                    av[ii] += sgn * dv[jf - 1] * fpow_int(b, mm);
                    at[jf - 1] = dv[jf - 1];
                }
                av[ii] = b * av[ii] * fpow_int(epai, n) / nf;
            }
        }

        // Build f(r) = r^3 / (1 + exp((r-r_nuc)/epai))
        // Fortran l tracks last i (1-based) where dv was set.
        // We use l_f to track the Fortran 1-based index.
        int l_f = 0;
        for (int i = 0; i < np; ++i) {
            b = 1.0 + std::exp((dr[i] - dr[nuc0]) / epai);
            if (b * av[3] > 1.0e+15) break;  // goto 51
            dv[i] = dr[i] * dr[i] * dr[i] / b;
            l_f = i + 1;  // 1-based Fortran index
        }
        // Fortran: if (l >= np-1) l = np-2
        if (l_f >= np - 1) l_f = np - 2;
        // Convert to 0-based for C++ array access
        int l0 = l_f - 1;  // 0-based index of last valid point
        int k = l_f;       // Fortran k = l+1 (1-based), which is l_f+1-1 = l_f in 0-based
        for (int i = k; i < np; ++i) {
            dv[i] = 0.0;
        }

        // Integration setup
        at[0] = 0.0;
        at[1] = 0.0;
        int kk = 2;
        for (int i = 3; i < ndor; ++i) {
            // Fortran: i goes 4..ndor, k increments from 3
            kk++;
            for (int jj = 0; jj < 2; ++jj) {
                at[jj] += av[i] * fpow_int(dr[jj], kk) / kk;
            }
            av[i] = av[i] / (kk * (kk - 1));
            av[1] += av[i] * fpow_int(dr[0], kk);
        }

        // 4-point integration
        a = hx / 24.0;
        b = a * 13.0;
        int k2 = l0 + 1;
        for (int i = 2; i <= k2; ++i) {
            at[i] = at[i - 1] + b * (dv[i - 1] + dv[i]) - a * (dv[i - 2] + dv[i + 1]);
        }
        dv[l0] = at[l0];
        for (int i = k2; i < np; ++i) {
            dv[i] = dv[l0];
        }

        double e = std::exp(hx);
        double c = 1.0 / (e * e);

        // Backward loop (Fortran goto 83 loop)
        for (int i = l0 - 1; i >= 1; --i) {
            dv[i] = dv[i + 1] / e + b * (at[i + 1] / e + at[i]) - a * (at[i + 2] * c + at[i - 1] * e);
        }

        // Fortran: dv(1) = dv(3)*c + hx*(at(1)+4*at(2)/e+at(3)*c)/3
        dv[0] = dv[2] * c + hx * (at[0] + 4.0 * at[1] / e + at[2] * c) / 3.0;

        // av(2) = (av(2)+dv(1))/dr(1)
        av[1] = (av[1] + dv[0]) / dr[0];

        // a = -dz/dv(l)  (Fortran l is the last valid index, 1-based)
        a = -dz / dv[l0];

        for (int i = 3; i < ndor; ++i) {
            av[i] = -a * av[i];
        }
        av[1] = a * av[1];
        for (int i = 0; i < np; ++i) {
            dv[i] = a * dv[i] / dr[i];
        }
    }
}

} // namespace feff::atom
