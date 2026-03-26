// Quantum Debye-Waller factor using correlated Debye model.
// Converted from: src/DEBYE/sigms.f
// Coded by J. Rehr (29 July 91), subroutine version Dec 1991 (S. Zabinsky).
//
// The Debye model correlation function is:
//   c(Ri,Rj) = (3kT / mu*omega_D^2) * sqrt(mu^2 / mi*mj) * I
// where I = integral_0^1 (y/x) dw sin(w*x) * coth(w*y/2)
// and x = k_D*R, y = thetaD/T.

#include "debye.hpp"

#include <feff/constants.hpp>

#include "../common/periodic_table.hpp"

#include "../math/distance.hpp"

#include <cmath>

namespace feff::debye {

// Thread-local state for the integrand (replaces Fortran COMMON /xy/)
static thread_local double g_x = 0.0;
static thread_local double g_yinv = 0.0;

// Quantum Debye integrand: fn(w) = sin(w*x)/x * coth(w*y/2)
// At w=0: fn = 2*yinv (L'Hopital limit)
static double fn(double w) {
    double result = 2.0 * g_yinv;
    if (w < 1.0e-20) return result;

    double fac = w;
    if (g_x > 0.0) fac = std::sin(w * g_x) / g_x;

    double emwy = std::exp(-w / g_yinv);
    result = fac * (1.0 + emwy) / (1.0 - emwy);
    return result;
}

void bingrt(double& b, double& eps, int& n) {
    // Romberg integration of fn over [0,1]
    // Error is approximately 2^(-2n) ~ 10^(-0.6n)
    constexpr int nmax = 10;
    constexpr double tol = 1.0e-5;

    n = 0;
    int itn = 1;
    double del = 1.0;
    double bn = (fn(0.0) + fn(1.0)) / 2.0;
    double bo = bn;

    for (;;) {
        n++;
        if (n > nmax) {
            // Not converged warning
            return;
        }
        del /= 2.0;
        double sum = 0.0;
        for (int i = 1; i <= itn; ++i) {
            double zi = (2 * i - 1) * del;
            sum += fn(zi);
        }

        // bnp1 = b_{n+1} current value
        double bnp1 = bn / 2.0 + del * sum;
        // Cancel leading error: b = (4*bnp1 - bn) / 3
        b = (4.0 * bnp1 - bn) / 3.0;
        eps = std::abs((b - bo) / b);
        if (eps < tol) return;

        bn = bnp1;
        bo = b;
        itn *= 2;
    }
}

// hbar^2 / (kB * amu) * 10^20 in Angstrom^2 units
// Updated 2017: hbar = 1.054571800e-34, amu = 1.660539040e-27, kB = 1.38064852e-23
static constexpr double con = 48.50875019927435;

void corrfn(double rij, double& cij, double thetad, double tk,
            int iz1, int iz2, double rsavg) {
    double ami = common::atomic_weights[iz1];
    double amj = common::atomic_weights[iz2];
    double rs = rsavg;

    g_yinv = tk / thetad;
    double xkd = std::pow(9.0 * pi / 2.0, third) / (rs * bohr);
    double fac = (3.0 / 2.0) * con / (thetad * std::sqrt(ami * amj));
    g_x = xkd * rij;

    double grater = 0.0, eps_out = 0.0;
    int nx = 0;
    bingrt(grater, eps_out, nx);
    cij = fac * grater;
}

void sigms(double tk, double thetad, double rs, int nlegx, int nleg,
           const double rat[][3], const int iz[], double& sig2) {

    double sigtot = 0.0;

    for (int il = 1; il <= nleg; ++il) {
        for (int jl = il; jl <= nleg; ++jl) {
            // Distances between leg endpoints
            double rij = math::dist(rat[il], rat[jl]);
            double cij = 0.0;
            corrfn(rij, cij, thetad, tk, iz[il], iz[jl], rs);
            double sig2ij = cij;

            double rimjm = math::dist(rat[il - 1], rat[jl - 1]);
            double cimjm = 0.0;
            corrfn(rimjm, cimjm, thetad, tk, iz[il - 1], iz[jl - 1], rs);
            sig2ij += cimjm;

            double rijm = math::dist(rat[il], rat[jl - 1]);
            double cijm = 0.0;
            corrfn(rijm, cijm, thetad, tk, iz[il], iz[jl - 1], rs);
            sig2ij -= cijm;

            double rimj = math::dist(rat[il - 1], rat[jl]);
            double cimj = 0.0;
            corrfn(rimj, cimj, thetad, tk, iz[il - 1], iz[jl], rs);
            sig2ij -= cimj;

            double riim = math::dist(rat[il], rat[il - 1]);
            double rjjm = math::dist(rat[jl], rat[jl - 1]);

            // Dot product of leg vectors
            double ridotj = 0.0;
            for (int k = 0; k < 3; ++k) {
                ridotj += (rat[il][k] - rat[il - 1][k]) *
                          (rat[jl][k] - rat[jl - 1][k]);
            }
            ridotj /= (riim * rjjm);

            // Double count i != j terms
            if (jl != il) sig2ij *= 2.0;
            sig2ij *= ridotj;
            sigtot += sig2ij;
        }
    }

    sig2 = sigtot / 4.0;
    // sig2 is in Angstrom^2
}

} // namespace feff::debye
