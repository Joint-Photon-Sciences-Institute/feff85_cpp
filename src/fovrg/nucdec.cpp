// Nuclear potential construction for FOVRG module.
// Converted from: src/FOVRG/nucdec.f
//
// Constructs the nuclear potential and its development coefficients at origin.
// Supports point-charge nucleus, uniform distribution, and Fermi distribution.

#include "nucdec.hpp"
#include <feff/dimensions.hpp>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace feff::fovrg {

static inline double fpow_int(double x, int n) {
    if (n == 0) return 1.0;
    if (n < 0) { x = 1.0/x; n = -n; }
    double r = 1.0;
    for (int i = 0; i < n; ++i) r *= x;
    return r;
}

void nucdec(double av[10], double dr[], double dv[], double dz,
            double hx, int& nuc, int np, int ndor, double& dr1)
{
    // Local workspace
    double at[nrptx];

    // Atomic mass and Fermi distribution parameter
    // a = 0 means point charge; epai = 0 means uniform distribution
    double a = 0.0;
    double epai = 0.0;

    // Calculate radial mesh
    if (a <= 1.0e-01) {
        nuc = 1;
    } else {
        // Fortran: a**(1./3.) — single-precision exponent (1./3. = 0.333333343267...)
        // Currently dead code (a=0 above), but match Fortran for correctness.
        a = dz * std::pow(a, static_cast<double>(1.0f / 3.0f)) * 2.2677e-05;
        double b = a / std::exp(hx * (nuc - 1));
        if (b <= dr1) {
            dr1 = b;
        } else {
            b = std::log(a / dr1) / hx;
            nuc = 3 + 2 * static_cast<int>(b / 2.0);
            if (nuc >= np) {
                throw std::runtime_error("dr1 too small");
            }
            dr1 = a * std::exp(-(nuc - 1) * hx);
        }
    }

    // Build radial mesh: dr[i] = dr1/dz * exp(hx * i) for i=0..np-1
    // Fortran: dr(1) = dr1/dz, dr(l) = dr(1)*exp(hx*(l-1))
    dr[0] = dr1 / dz;
    for (int l = 1; l < np; l++) {
        dr[l] = dr[0] * std::exp(hx * l);
    }

    if (ndor < 5) {
        std::cerr << "stopped in nucdec, ndor should be > 4." << std::endl;
        throw std::runtime_error("NUCDEC-1");
    }

    // Zero development coefficients
    for (int i = 0; i < ndor; i++) {
        av[i] = 0.0;
    }

    if (epai <= 0.0) {
        // Uniform or point-charge distribution
        // dv[i] = -dz / dr[i]
        for (int i = 0; i < np; i++) {
            dv[i] = -dz / dr[i];
        }
        if (nuc <= 1) {
            // Point charge: av(1) = -dz
            av[0] = -dz;
        } else {
            // Uniform distribution inside nuclear radius
            // Fortran: av(2), av(4) with 1-based indexing
            av[1] = -3.0 * dz / (dr[nuc - 1] + dr[nuc - 1]);
            av[3] = -av[1] / (3.0 * dr[nuc - 1] * dr[nuc - 1]);
            int l = nuc - 1;  // Fortran nuc-1, but 0-based means indices 0..nuc-2
            for (int i = 0; i < l; i++) {
                dv[i] = av[1] + av[3] * dr[i] * dr[i];
            }
        }
    } else {
        // Fermi distribution
        double b = std::exp(-dr[nuc - 1] / epai);
        b = 1.0 / (1.0 + b);
        av[3] = b;       // av(4) in Fortran
        av[4] = epai * b * (b - 1.0);  // av(5) in Fortran
        if (ndor > 5) {
            at[0] = 1.0;
            at[1] = 1.0;
            int nf = 1;
            for (int i = 5; i < ndor; i++) {
                // i is 0-based C++ index for av(i+1) in Fortran, which is av(6)..av(ndor)
                // Fortran loop: do i=6,ndor => C++ i=5..ndor-1
                int n = i - 3;  // Fortran: n = i-4, but i is 1-based there => n = (i+1)-4 = i-3
                nf = n * nf;
                dv[0] = n * at[0];
                int n1 = n + 1;
                dv[n1 - 1] = 1.0;  // dv(n1) in Fortran, 0-based: dv[n1-1]
                for (int j = 1; j < n; j++) {
                    // Fortran: do j=2,n => C++ j=1..n-1
                    dv[j] = (n - j + 1) * at[j - 1] + (n - j) * at[j];
                }
                for (int j = 0; j < n1; j++) {
                    // Fortran: do j=1,n1
                    int m = n - j;  // Fortran: m=n+1-j, but j is 1-based => m = n+1-(j+1) = n-j
                    int l = 1;
                    if ((j + 1) % 2 == 0) l = -l;  // Fortran: mod(j,2) with 1-based j
                    av[i] = av[i] + l * dv[j] * fpow_int(b, m);
                    at[j] = dv[j];
                }
                av[i] = b * av[i] * fpow_int(epai, n) / nf;
            }
        }

        int l = 0;
        for (int i = 0; i < np; i++) {
            b = 1.0 + std::exp((dr[i] - dr[nuc - 1]) / epai);
            if ((b * av[3]) > 1.0e+15) break;
            dv[i] = dr[i] * dr[i] * dr[i] / b;
            l = i;
        }
        // Fortran: if (l.ge.(np-1)) l=np-2
        if (l >= np - 1) l = np - 2;
        int k = l + 1;
        for (int i = k; i < np; i++) {
            dv[i] = 0.0;
        }
        at[0] = 0.0;
        at[1] = 0.0;
        k = 2;
        for (int i = 3; i < ndor; i++) {
            // Fortran: do i=4,ndor => 0-based i=3..ndor-1
            k = k + 1;
            for (int j = 0; j < 2; j++) {
                at[j] = at[j] + av[i] * fpow_int(dr[j], k) / k;
            }
            av[i] = av[i] / (k * (k - 1));
            av[1] = av[1] + av[i] * fpow_int(dr[0], k);
        }
        a = hx / 24.0;
        b = a * 13.0;
        k = l + 1;
        for (int i = 2; i <= k; i++) {
            // Fortran: do i=3,k => 0-based i=2..k
            // at(i) = at(i-1) + b*(dv(i-1)+dv(i)) - a*(dv(i-2)+dv(i+1))
            at[i] = at[i - 1] + b * (dv[i - 1] + dv[i]) - a * (dv[i - 2] + dv[i + 1]);
        }
        dv[l] = at[l];
        for (int i = k; i < np; i++) {
            dv[i] = dv[l];
        }
        double e = std::exp(hx);
        double c = 1.0 / (e * e);
        int ii = l - 1;
        while (ii > 0) {
            // Fortran: dv(i)=dv(i+1)/e+b*(at(i+1)/e+at(i))-a*(at(i+2)*c+at(i-1)*e)
            dv[ii] = dv[ii + 1] / e + b * (at[ii + 1] / e + at[ii]) - a * (at[ii + 2] * c + at[ii - 1] * e);
            ii = ii - 1;
        }
        // Fortran: dv(1) = dv(3)*c + hx*(at(1)+4*at(2)/e+at(3)*c)/3
        dv[0] = dv[2] * c + hx * (at[0] + 4.0 * at[1] / e + at[2] * c) / 3.0;
        av[1] = (av[1] + dv[0]) / dr[0];
        a = -dz / dv[l];
        for (int i = 3; i < ndor; i++) {
            av[i] = -a * av[i];
        }
        av[1] = a * av[1];
        for (int i = 0; i < np; i++) {
            dv[i] = a * dv[i] / dr[i];
        }
    }
}

} // namespace feff::fovrg
