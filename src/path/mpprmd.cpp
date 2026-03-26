// Double-precision path parameter computation for output.
// Converted from: src/PATH/mpprmd.f
// Contains strigd and sargd helper functions as local statics.

#include "mpprmd.hpp"
#include "../math/distance.hpp"
#include <cmath>
#include <complex>

namespace feff::path {

// Helper: compute cos(theta), sin(theta), cos(phi), sin(phi) from (x,y,z)
static void strigd(double x, double y, double z,
                   double& ct, double& st, double& cp, double& sp) {
    constexpr double eps = 1.0e-6;
    double r = std::sqrt(x*x + y*y + z*z);
    double rxy = std::sqrt(x*x + y*y);

    if (r < eps) {
        ct = 1.0; st = 0.0;
    } else {
        ct = z / r;
        st = rxy / r;
    }
    if (rxy < eps) {
        cp = 1.0; sp = 0.0;
    } else {
        cp = x / rxy;
        sp = y / rxy;
    }
}

// Helper: extract angle from complex exponential
static void sargd(std::complex<double> c, double& th) {
    constexpr double eps = 1.0e-6;
    double x = c.real();
    double y = c.imag();
    if (std::abs(x) < eps) x = 0.0;
    if (std::abs(y) < eps) y = 0.0;
    if (std::abs(x) < eps && std::abs(y) < eps) {
        th = 0.0;
    } else {
        th = std::atan2(y, x);
    }
}

void mpprmd(int npat, const int ipat[], double ri[], double beta[], double eta[],
            const AtomData& atoms) {

    const std::complex<double> coni_local(0.0, 1.0);

    int n = npat + 1;  // nleg

    std::vector<std::complex<double>> alph(n);
    std::vector<std::complex<double>> gamm(n + 1); // need gamm[n] = gamm(npat+2) = gamm(1)

    for (int jj = 0; jj < n; ++jj) {
        int j = jj;  // 0-based

        // Determine atoms i, ip1, im1
        int i_atom, ip1, im1;
        if (j == n - 1) {
            // j is central atom (last leg)
            i_atom = 0;
            ip1 = ipat[0];
            im1 = ipat[npat - 1];
        } else if (j == npat - 1) {
            // j is last path atom
            i_atom = ipat[j];
            ip1 = 0;
            im1 = (npat == 1) ? 0 : ipat[npat - 2];
        } else if (j == 0) {
            // j is first atom
            i_atom = ipat[j];
            ip1 = (npat == 1) ? 0 : ipat[j + 1];
            im1 = 0;
        } else {
            i_atom = ipat[j];
            ip1 = ipat[j + 1];
            im1 = ipat[j - 1];
        }

        // Compute spherical angles for forward direction (i -> ip1)
        double x = atoms.rat[ip1][0] - atoms.rat[i_atom][0];
        double y = atoms.rat[ip1][1] - atoms.rat[i_atom][1];
        double z = atoms.rat[ip1][2] - atoms.rat[i_atom][2];
        double ct, st, cp, sp;
        strigd(x, y, z, ct, st, cp, sp);

        // Compute spherical angles for backward direction (im1 -> i)
        x = atoms.rat[i_atom][0] - atoms.rat[im1][0];
        y = atoms.rat[i_atom][1] - atoms.rat[im1][1];
        z = atoms.rat[i_atom][2] - atoms.rat[im1][2];
        double ctp, stp, cpp, spp;
        strigd(x, y, z, ctp, stp, cpp, spp);

        // cppp = cos(phi' - phi), sppp = sin(phi' - phi)
        double cppp = cp * cpp + sp * spp;
        double sppp = spp * cp - cpp * sp;

        // Euler angle components
        alph[jj] = st * ctp - ct * stp * cppp - coni_local * stp * sppp;
        beta[jj] = ct * ctp + st * stp * cppp;
        // Clamp beta to [-1, 1]
        if (beta[jj] < -1.0) beta[jj] = -1.0;
        if (beta[jj] >  1.0) beta[jj] =  1.0;
        gamm[jj] = st * ctp * cppp - ct * stp + coni_local * st * sppp;

        // Leg distance
        float r0[3] = {atoms.rat[i_atom][0], atoms.rat[i_atom][1], atoms.rat[i_atom][2]};
        float r1[3] = {atoms.rat[im1][0], atoms.rat[im1][1], atoms.rat[im1][2]};
        ri[jj] = feff::math::sdist(r0, r1);
    }

    // eta(j) = arg(alph(j) * gamm(j+1)), with gamm(npat+2) = gamm(1)
    gamm[n] = gamm[0];
    for (int j = 0; j < n; ++j) {
        std::complex<double> eieta = alph[j] * gamm[j + 1];
        sargd(eieta, eta[j]);
    }

    // Convert beta from cos(beta) to angle
    for (int j = 0; j < n; ++j) {
        if (beta[j] >  1.0) beta[j] =  1.0;
        if (beta[j] < -1.0) beta[j] = -1.0;
        beta[j] = std::acos(beta[j]);
    }
}

} // namespace feff::path
