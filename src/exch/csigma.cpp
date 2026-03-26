// Many-pole self-energy with self-consistent momentum
// Converted from src/EXCH/csigma.f
//
// Implements the Hedin-Lundqvist self-energy model including:
//   - Self-consistent momentum iteration
//   - Adaptive Gaussian quadrature (cgratr)
//   - Integration kernels (r1, r2, r3, dr1, dr2, dr3)
//   - Plasmon dispersion relation (fq)
//   - Hartree-Fock exchange (hfexc)

#include "csigma.hpp"
#include "fndsng.hpp"

#include <cmath>
#include <algorithm>
#include <stdexcept>

#include <feff/types.hpp>
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>

namespace feff::exch {

// ============================================================================
// Plasmon dispersion relation
// ============================================================================

FeffComplex fq_dispersion(FeffComplex q, const double dppar[10])
{
    double xwg  = dppar[0];
    double xgam = dppar[1];

    // fq(q) = sqrt((wp - i*gamma)^2 + 4/3*q^2 + q^4)
    // electron gas parameters: a2=4/3, a4=1
    constexpr double a2 = 4.0 / 3.0;
    constexpr double a4 = 1.0;

    FeffComplex wp_c(xwg, -xgam);
    FeffComplex result = wp_c * wp_c + a2 * q * q + a4 * q * q * q * q;
    return std::sqrt(result);
}

// ============================================================================
// Integration kernels
// ============================================================================

FeffComplex kernel_r1(FeffComplex q, const double dppar[10], const FeffComplex cpar[10])
{
    FeffComplex ck  = cpar[0];
    FeffComplex xe  = cpar[1];
    double      xeg = dppar[3];

    FeffComplex fqq = fq_dispersion(q, dppar);
    FeffComplex fiq = 1.0 / (q * fqq);

    FeffComplex a1 = 1.0 - xeg - xe - fqq - coni * 1.0e-10;
    FeffComplex a2 = (ck + q) * (ck + q) - xe + fqq - coni * 1.0e-10;
    FeffComplex a3 = (ck - q) * (ck - q) - xe - fqq - coni * 1.0e-10;
    FeffComplex a4 = 1.0 + xeg - xe + fqq - coni * 1.0e-10;

    return fiq * (std::log(a1) + std::log(a2) - std::log(a3) - std::log(a4));
}

FeffComplex kernel_r2(FeffComplex q, const double dppar[10], const FeffComplex cpar[10])
{
    FeffComplex ck = cpar[0];
    FeffComplex xe = cpar[1];

    FeffComplex fqq = fq_dispersion(q, dppar);
    FeffComplex fiq = 1.0 / (q * fqq);

    FeffComplex a1 = (ck + q) * (ck + q) - xe + fqq - coni * 1.0e-10;
    FeffComplex a2 = (ck - q) * (ck - q) - xe + fqq - coni * 1.0e-10;

    return fiq * (std::log(a1) - std::log(a2));
}

FeffComplex kernel_r3(FeffComplex q, const double dppar[10], const FeffComplex cpar[10])
{
    FeffComplex ck = cpar[0];
    FeffComplex xe = cpar[1];

    // valid only for k < kF, q < kF - k
    FeffComplex fqq = fq_dispersion(q, dppar);
    FeffComplex fiq = 1.0 / (q * fqq);

    FeffComplex a1 = (ck + q) * (ck + q) - xe - fqq - coni * 1.0e-10;
    FeffComplex a2 = (ck - q) * (ck - q) - xe - fqq - coni * 1.0e-10;

    return fiq * (std::log(a1) - std::log(a2));
}

FeffComplex kernel_dr1(FeffComplex q, const double dppar[10], const FeffComplex cpar[10])
{
    FeffComplex ck  = cpar[0];
    FeffComplex xe  = cpar[1];
    double      xeg = dppar[3];

    FeffComplex fqq = fq_dispersion(q, dppar);
    FeffComplex fiq = 1.0 / (q * fqq);

    FeffComplex a1 = 1.0 - xeg - xe - fqq - coni * 1.0e-10;
    FeffComplex a2 = (ck + q) * (ck + q) - xe + fqq - coni * 1.0e-10;
    FeffComplex a3 = (ck - q) * (ck - q) - xe - fqq - coni * 1.0e-10;
    FeffComplex a4 = 1.0 + xeg - xe + fqq - coni * 1.0e-10;

    return -fiq * (1.0 / a1 + 1.0 / a2 - 1.0 / a3 - 1.0 / a4);
}

FeffComplex kernel_dr2(FeffComplex q, const double dppar[10], const FeffComplex cpar[10])
{
    FeffComplex ck = cpar[0];
    FeffComplex xe = cpar[1];

    FeffComplex fqq = fq_dispersion(q, dppar);
    FeffComplex fiq = 1.0 / (q * fqq);

    FeffComplex a1 = (ck + q) * (ck + q) - xe + fqq - coni * 1.0e-10;
    FeffComplex a2 = (ck - q) * (ck - q) - xe + fqq - coni * 1.0e-10;

    return -fiq * (1.0 / a1 - 1.0 / a2);
}

FeffComplex kernel_dr3(FeffComplex q, const double dppar[10], const FeffComplex cpar[10])
{
    FeffComplex ck = cpar[0];
    FeffComplex xe = cpar[1];

    // valid only for k < kF, q < kF - k
    FeffComplex fqq = fq_dispersion(q, dppar);
    FeffComplex fiq = 1.0 / (q * fqq);

    FeffComplex a1 = (ck + q) * (ck + q) - xe - fqq - coni * 1.0e-10;
    FeffComplex a2 = (ck - q) * (ck - q) - xe - fqq - coni * 1.0e-10;

    return -fiq * (1.0 / a1 - 1.0 / a2);
}

// ============================================================================
// Adaptive Gaussian quadrature (cgratr)
// ============================================================================

FeffComplex cgratr(KernelFn fn, const double dppar[10], const FeffComplex cpar[10],
                   FeffComplex xmin, FeffComplex xmax,
                   double abr, double rlr,
                   int nsing, FeffComplex xsing[20],
                   double& error, int& numcal, int& maxns)
{
    constexpr int mx = 1500;

    // Gauss-Legendre 3-point nodes and weights
    static constexpr double dx[3] = {
        0.1127016653792583,
        0.5,
        0.8872983346207417
    };
    static constexpr double wt[3] = {
        0.277777777777777778,
        0.4444444444444444444,
        0.2777777777777777778
    };
    // 9-point composite weights (3 x 3-point combined)
    static constexpr double wt9[9] = {
        0.0616938806304841571,
        0.108384229110206161,
        0.0398463603260281088,
        0.175209035316976464,
        0.229732989232610220,
        0.175209035316976464,
        0.0398463603260281088,
        0.108384229110206161,
        0.0616938806304841571
    };

    FeffComplex fval[3][mx];
    FeffComplex xleft[mx];

    // nstack is the number of different intervals
    int nstack = nsing + 1;
    maxns = nstack;
    error = 0.0;
    FeffComplex result(0.0, 0.0);

    // Initial placement of regions bounded by singularities
    xleft[0] = xmin;
    xleft[nsing + 1] = xmax;
    for (int j = 0; j < nsing; ++j) {
        xleft[j + 1] = xsing[j];
    }

    // For each region, calculate the function at three selected points
    for (int jj = 0; jj < nstack; ++jj) {
        FeffComplex del = xleft[jj + 1] - xleft[jj];
        for (int j = 0; j < 3; ++j) {
            fval[j][jj] = fn(xleft[jj] + del * dx[j], dppar, cpar);
        }
    }
    numcal = nstack * 3;

    // Main adaptive loop
    for (;;) {
        if (nstack + 3 >= mx) {
            throw std::runtime_error("cgratr: too many regions");
        }

        // Divide the rightmost region into three subregions
        FeffComplex del = xleft[nstack] - xleft[nstack - 1];
        xleft[nstack + 2] = xleft[nstack];
        xleft[nstack]     = xleft[nstack - 1] + del * dx[0] * 2.0;
        xleft[nstack + 1] = xleft[nstack + 2] - del * dx[0] * 2.0;

        // Redistribute the three existing data points
        fval[1][nstack + 1] = fval[2][nstack - 1];
        fval[1][nstack]     = fval[1][nstack - 1];
        fval[1][nstack - 1] = fval[0][nstack - 1];

        // Compute value (9-point) and valu (3x3-point) over the three subregions
        int icount = 0;
        FeffComplex value(0.0, 0.0);
        FeffComplex valu(0.0, 0.0);

        for (int j = nstack - 1; j <= nstack + 1; ++j) {
            FeffComplex del1 = xleft[j + 1] - xleft[j];
            fval[0][j] = fn(xleft[j] + dx[0] * del1, dppar, cpar);
            fval[2][j] = fn(xleft[j] + dx[2] * del1, dppar, cpar);
            numcal += 2;

            for (int k = 0; k < 3; ++k) {
                value += wt9[icount] * fval[k][j] * del;
                valu  += fval[k][j] * wt[k] * del1;
                ++icount;
            }
        }

        double dif = std::abs(value - valu);
        double frac = std::abs(del / (xmax - xmin));
        bool atsing = (frac <= 1.0e-8);

        if (dif <= abr * frac || dif <= rlr * std::abs(value) ||
            (atsing && (frac <= 1.0e-15 || dif <= abr * 0.1))) {
            // Accurate enough: add to total and pop the rightmost region
            result += value;
            error += std::abs(dif);
            --nstack;
            if (nstack <= 0) return result;
        } else {
            // Not accurate enough: keep the three subregions
            nstack += 2;
            maxns = std::max(maxns, nstack);
        }
    }

    return result; // unreachable
}

// ============================================================================
// Hartree-Fock exchange
// ============================================================================

FeffComplex hfexc(FeffComplex ck_in, double efermi, double kfermi)
{
    FeffComplex ck = ck_in / kfermi;
    FeffComplex c = -2.0 * efermi / (pi * kfermi);

    if (std::abs(ck - 1.0) <= 0.00001) {
        return c;
    } else {
        return c * (1.0 + (1.0 / ck - ck) * std::log((1.0 + ck) / (ck - 1.0)) / 2.0);
    }
}

// ============================================================================
// Sigma1: single-pole energy-dependent self-energy
// ============================================================================

FeffComplex sigma1(FeffComplex ck, FeffComplex energy,
                   double wi, double gamma, double amp,
                   double kfermi, double efermi)
{
    constexpr double zero_pl = 1.0e-5;
    constexpr double inf     = 1.0e2;
    constexpr double abs_err = 1.0e-5;
    constexpr double rel_err = 1.0e-4;

    int nsing = 0, ncalls = 0, maxnr = 0;
    double dppar[10] = {};
    FeffComplex cpar[10] = {};
    FeffComplex xsing[20] = {};
    double error;

    // Set up dimensionless parameters
    dppar[0] = wi / efermi;           // xwg
    dppar[1] = gamma / efermi;        // xgam
    dppar[2] = energy.real() / efermi; // xe
    dppar[3] = 0.0;                   // xeg (gap energy)

    cpar[0] = ck / kfermi;            // dimensionless ck
    cpar[1] = energy / efermi;        // complex energy

    FeffComplex result(0.0, 0.0);

    // Loop i1 = 1 to 1 (only one iteration in current code)
    {
        // 1) Integral of r2 from ck/kF + 1 to Inf
        FeffComplex limit1 = ck / kfermi + 1.0;
        FeffComplex limit2(inf, 0.0);

        int ifcn = 2;
        fndsng(limit1, limit2, nsing, xsing, dppar, cpar, ifcn);

        FeffComplex hlint1 = cgratr(kernel_r2, dppar, cpar, limit1, limit2,
                                     abs_err, rel_err, nsing, xsing,
                                     error, ncalls, maxnr);
        // Reset xsing imaginary parts
        for (int i = 0; i < nsing; ++i) {
            xsing[i] = FeffComplex(xsing[i].real(), 0.0);
        }

        // 2) Integral of r1 from |ck/kF - 1| to ck/kF + 1
        limit1 = FeffComplex(std::max(std::abs(ck.real() / kfermi - 1.0), zero_pl), 0.0);
        limit2 = ck / kfermi + 1.0;

        ifcn = 1;
        fndsng(limit1, limit2, nsing, xsing, dppar, cpar, ifcn);

        FeffComplex hlint2 = cgratr(kernel_r1, dppar, cpar, limit1, limit2,
                                     abs_err, rel_err, nsing, xsing,
                                     error, ncalls, maxnr);
        for (int i = 0; i < nsing; ++i) {
            xsing[i] = FeffComplex(xsing[i].real(), 0.0);
        }

        // 3) Third integral depending on ck vs kF
        FeffComplex hlint3(0.0, 0.0);
        limit1 = FeffComplex(zero_pl, 0.0);
        limit2 = FeffComplex(std::abs(ck.real() / kfermi - 1.0), 0.0);

        if (std::abs(ck.real() - kfermi) < zero_pl ||
            limit2.real() <= limit1.real()) {
            // ck == kF: hlint3 = 0
            hlint3 = FeffComplex(0.0, 0.0);
        } else if (ck.real() < kfermi) {
            // ck < kF: integrate r3 from 0 to 1 - ck/kF
            limit2 = FeffComplex(1.0 - ck.real() / kfermi, 0.0);
            ifcn = 3;
            fndsng(limit1, limit2, nsing, xsing, dppar, cpar, ifcn);
            hlint3 = cgratr(kernel_r3, dppar, cpar, limit1, limit2,
                            abs_err, rel_err, nsing, xsing,
                            error, ncalls, maxnr);
            for (int i = 0; i < nsing; ++i) {
                xsing[i] = FeffComplex(xsing[i].real(), 0.0);
            }
        } else {
            // ck > kF: integrate r2 from 0 to ck/kF - 1
            limit2 = FeffComplex(ck.real() / kfermi - 1.0, 0.0);
            ifcn = 2;
            fndsng(limit1, limit2, nsing, xsing, dppar, cpar, ifcn);
            hlint3 = cgratr(kernel_r2, dppar, cpar, limit1, limit2,
                            abs_err, rel_err, nsing, xsing,
                            error, ncalls, maxnr);
            for (int i = 0; i < nsing; ++i) {
                xsing[i] = FeffComplex(xsing[i].real(), 0.0);
            }
        }

        // Sigma1 = -Amp * Wi * (Wi - i*Gamma) / (2*pi*EFermi*ck) * (sum of integrals)
        result = -amp * wi * (wi - coni * gamma) / (2.0 * pi * efermi * ck)
                 * (hlint1 + hlint2 + hlint3);
    }

    return result;
}

// ============================================================================
// dSigma: derivative of self-energy w.r.t. energy
// ============================================================================

FeffComplex dsigma(FeffComplex ck, FeffComplex energy,
                   double wi, double gamma, double amp,
                   double kfermi, double efermi)
{
    constexpr double zero_pl = 1.0e-5;
    constexpr double inf     = 1.0e2;
    constexpr double abs_err = 1.0e-5;
    constexpr double rel_err = 1.0e-4;

    int nsing = 0, ncalls = 0, maxnr = 0;
    double dppar[10] = {};
    FeffComplex cpar[10] = {};
    FeffComplex xsing[20] = {};
    double error;

    // Set up dimensionless parameters
    dppar[0] = wi / efermi;
    dppar[1] = gamma / efermi;
    dppar[2] = energy.real() / efermi;
    dppar[3] = 0.0;

    cpar[0] = ck / kfermi;
    // For derivative, cpar[1] includes broadening: E/EF + i*Gamma/EF
    cpar[1] = energy / efermi + coni * dppar[1];

    FeffComplex result(0.0, 0.0);

    // Single iteration (i1 = 1)
    {
        // 1) Integral of dr2 from ck/kF + 1 to Inf
        FeffComplex limit1 = ck / kfermi + 1.0;
        FeffComplex limit2(inf, 0.0);

        int ifcn = 2;
        fndsng(limit1, limit2, nsing, xsing, dppar, cpar, ifcn);

        FeffComplex hlint1 = cgratr(kernel_dr2, dppar, cpar, limit1, limit2,
                                     abs_err, rel_err, nsing, xsing,
                                     error, ncalls, maxnr);
        for (int i = 0; i < nsing; ++i) {
            xsing[i] = FeffComplex(xsing[i].real(), 0.0);
        }

        // 2) Integral of dr1 from |ck/kF - 1| to ck/kF + 1
        limit1 = FeffComplex(std::max(std::abs(ck.real() / kfermi - 1.0), zero_pl), 0.0);
        limit2 = ck / kfermi + 1.0;

        ifcn = 1;
        fndsng(limit1, limit2, nsing, xsing, dppar, cpar, ifcn);

        FeffComplex hlint2 = cgratr(kernel_dr1, dppar, cpar, limit1, limit2,
                                     abs_err, rel_err, nsing, xsing,
                                     error, ncalls, maxnr);
        for (int i = 0; i < nsing; ++i) {
            xsing[i] = FeffComplex(xsing[i].real(), 0.0);
        }

        // 3) Third integral depending on ck vs kF
        FeffComplex hlint3(0.0, 0.0);
        limit1 = FeffComplex(zero_pl, 0.0);
        limit2 = FeffComplex(std::abs(ck.real() / kfermi - 1.0), 0.0);

        if (std::abs(ck.real() - kfermi) < zero_pl ||
            limit2.real() <= limit1.real()) {
            hlint3 = FeffComplex(0.0, 0.0);
        } else if (ck.real() < kfermi) {
            limit2 = FeffComplex(1.0 - ck.real() / kfermi, 0.0);
            ifcn = 3;
            fndsng(limit1, limit2, nsing, xsing, dppar, cpar, ifcn);
            hlint3 = cgratr(kernel_dr3, dppar, cpar, limit1, limit2,
                            abs_err, rel_err, nsing, xsing,
                            error, ncalls, maxnr);
            for (int i = 0; i < nsing; ++i) {
                xsing[i] = FeffComplex(xsing[i].real(), 0.0);
            }
        } else {
            limit2 = FeffComplex(ck.real() / kfermi - 1.0, 0.0);
            ifcn = 2;
            fndsng(limit1, limit2, nsing, xsing, dppar, cpar, ifcn);
            hlint3 = cgratr(kernel_dr2, dppar, cpar, limit1, limit2,
                            abs_err, rel_err, nsing, xsing,
                            error, ncalls, maxnr);
            for (int i = 0; i < nsing; ++i) {
                xsing[i] = FeffComplex(xsing[i].real(), 0.0);
            }
        }

        result = -amp * wi * (wi - coni * gamma) / (2.0 * pi * efermi * ck)
                 * (hlint1 + hlint2 + hlint3);
    }

    return result;
}

// ============================================================================
// CSigma: many-pole self-energy with self-consistent momentum
// ============================================================================

void csigma(FeffComplex energy, double mu, double rs,
            double& resig, double& imsig,
            const double wpscl[], const double ampfac[])
{
    constexpr int mx_iter = 1;

    double kfermi = fa / rs;
    double efermi = kfermi * kfermi / 2.0;

    FeffComplex sig_tot(0.0, 0.0);
    FeffComplex sigma_f(0.0, 0.0);
    double gam = 0.0;
    FeffComplex ckf(0.0, 0.0);
    FeffComplex ck0(0.0, 0.0);

    // Self-consistency loop (MxIter = 1, no actual iteration)
    for (int i2 = 0; i2 < mx_iter; ++i2) {

        // Loop over poles to find SigmaF (self-energy at Fermi level)
        for (int i1 = 0; i1 < MxPole; ++i1) {
            if (wpscl[i1] < -1000.0) break;

            double wp = std::sqrt(3.0 / (rs * rs * rs)) * wpscl[i1];

            ckf = FeffComplex(kfermi * 1.00001, 0.0);
            FeffComplex rel_en(efermi, 0.0);
            sigma_f += sigma1(ckf, rel_en, wp, gam, ampfac[i1], kfermi, efermi);
        }

        // Loop over poles to find Sigma0
        for (int i1 = 0; i1 < MxPole; ++i1) {
            if (wpscl[i1] < -1000.0) break;

            double wp = std::sqrt(3.0 / (rs * rs * rs)) * wpscl[i1];

            FeffComplex rel_en(energy.real() - mu + efermi, 0.0);
            ck0 = std::sqrt(2.0 * rel_en.real());

            FeffComplex sigma0 = sigma1(ck0, rel_en, wp, gam, ampfac[i1],
                                        kfermi, efermi);
            sig_tot += sigma0;
        }
    }

    // Form delta sigma
    sig_tot = sig_tot - sigma_f;

    // Add Hartree-Fock part
    FeffComplex del_hf = hfexc(ck0, efermi, kfermi) - hfexc(ckf, efermi, kfermi);
    sig_tot += del_hf;

    resig = sig_tot.real();
    imsig = sig_tot.imag();
}

} // namespace feff::exch
