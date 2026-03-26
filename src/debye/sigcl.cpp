// Classical Debye-Waller factor using classical Debye model.
// Converted from: src/DEBYE/sigcl.f
// Hacked for classical Debye Model, 7/2006, Kevin Jorissen.
//
// Same structure as sigms but uses the classical limit:
//   coth(w*y/2) -> 2/(w*y) = 2*yinv/w

#include "debye.hpp"

#include <feff/constants.hpp>

#include "../common/periodic_table.hpp"

#include "../math/distance.hpp"

#include <cmath>

namespace feff::debye {

// Thread-local state for classical integrand
static thread_local double g_x_cl = 0.0;
static thread_local double g_yinv_cl = 0.0;

// Classical Debye integrand: fn2(w) = sin(w*x)/x * 2*yinv/w
// At w=0: fn2 = 2*yinv (same limit)
static double fn2(double w) {
    double result = 2.0 * g_yinv_cl;
    if (w < 1.0e-20) return result;

    double fac = w;
    if (g_x_cl > 0.0) fac = std::sin(w * g_x_cl) / g_x_cl;

    // Classical limit: coth(w*y/2) -> 2/(w*y)
    result = fac * 2.0 * g_yinv_cl / w;
    return result;
}

void bingrt2(double& b, double& eps, int& n) {
    // Romberg integration of fn2 over [0,1]
    constexpr int nmax = 10;
    constexpr double tol = 1.0e-5;

    n = 0;
    int itn = 1;
    double del = 1.0;
    double bn = (fn2(0.0) + fn2(1.0)) / 2.0;
    double bo = bn;

    for (;;) {
        n++;
        if (n > nmax) return;

        del /= 2.0;
        double sum = 0.0;
        for (int i = 1; i <= itn; ++i) {
            double zi = (2 * i - 1) * del;
            sum += fn2(zi);
        }

        double bnp1 = bn / 2.0 + del * sum;
        b = (4.0 * bnp1 - bn) / 3.0;
        eps = std::abs((b - bo) / b);
        if (eps < tol) return;

        bn = bnp1;
        bo = b;
        itn *= 2;
    }
}

static constexpr double con = 48.50875019927435;

void corrfn2(double rij, double& cij, double thetad, double tk,
             int iz1, int iz2, double rsavg) {
    double ami = common::atomic_weights[iz1];
    double amj = common::atomic_weights[iz2];
    double rs = rsavg;

    g_yinv_cl = tk / thetad;
    double xkd = std::pow(9.0 * pi / 2.0, third) / (rs * bohr);
    double fac = (3.0 / 2.0) * con / (thetad * std::sqrt(ami * amj));
    g_x_cl = xkd * rij;

    double grater = 0.0, eps_out = 0.0;
    int nx = 0;
    bingrt2(grater, eps_out, nx);
    cij = fac * grater;
}

void sigcl(double tk, double thetad, double rs, int nlegx, int nleg,
           const double rat[][3], const int iz[], double& sig2) {

    double sigtot = 0.0;

    for (int il = 1; il <= nleg; ++il) {
        for (int jl = il; jl <= nleg; ++jl) {
            double rij = math::dist(rat[il], rat[jl]);
            double cij = 0.0;
            corrfn2(rij, cij, thetad, tk, iz[il], iz[jl], rs);
            double sig2ij = cij;

            double rimjm = math::dist(rat[il - 1], rat[jl - 1]);
            double cimjm = 0.0;
            corrfn2(rimjm, cimjm, thetad, tk, iz[il - 1], iz[jl - 1], rs);
            sig2ij += cimjm;

            double rijm = math::dist(rat[il], rat[jl - 1]);
            double cijm = 0.0;
            corrfn2(rijm, cijm, thetad, tk, iz[il], iz[jl - 1], rs);
            sig2ij -= cijm;

            double rimj = math::dist(rat[il - 1], rat[jl]);
            double cimj = 0.0;
            corrfn2(rimj, cimj, thetad, tk, iz[il - 1], iz[jl], rs);
            sig2ij -= cimj;

            double riim = math::dist(rat[il], rat[il - 1]);
            double rjjm = math::dist(rat[jl], rat[jl - 1]);

            double ridotj = 0.0;
            for (int k = 0; k < 3; ++k) {
                ridotj += (rat[il][k] - rat[il - 1][k]) *
                          (rat[jl][k] - rat[jl - 1][k]);
            }
            ridotj /= (riim * rjjm);

            if (jl != il) sig2ij *= 2.0;
            sig2ij *= ridotj;
            sigtot += sig2ij;
        }
    }

    sig2 = sigtot / 4.0;
}

} // namespace feff::debye
