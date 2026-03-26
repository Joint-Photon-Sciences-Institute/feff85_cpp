#pragma once

// Many-pole self-energy with self-consistent momentum
// Converted from src/EXCH/csigma.f (~894 lines)
//
// Calculates Sigma(k(Energy), Energy) based on an electron gas model
// of epsilon^{-1}, using the Hedin-Lundqvist formalism.
//
// Contains internal functions:
//   sigma1  - energy-dependent self-energy via adaptive Gaussian quadrature
//   dsigma  - derivative of self-energy w.r.t. energy
//   hfexc   - Hartree-Fock exchange (Dirac-Hara)
//   cgratr  - adaptive Gaussian quadrature for complex integrands
//   r1,r2,r3,dr1,dr2,dr3 - integration kernels
//   fq      - plasmon dispersion relation

#include <feff/types.hpp>

namespace feff::exch {

/// Many-pole self-energy with self-consistent momentum.
///
/// Solves: k0^2 = 2*Energy - 2*Mu - 2*(Sigma(k0,Energy)-Sigma(kFermi,EFermi))
///
/// @param energy  Energy at which to evaluate Sigma (complex, a.u.)
/// @param mu      Fermi energy from FEFF (a.u.)
/// @param rs      Density parameter (sphere of radius Rs holds charge e)
/// @param resig   [out] Re[Sigma(Energy, k(Energy))]
/// @param imsig   [out] Im[Sigma(Energy, k(Energy))]
/// @param wpscl   Plasmon pole frequency scale factors (MxPole array)
/// @param ampfac  Plasmon pole amplitude factors (MxPole array)
void csigma(FeffComplex energy, double mu, double rs,
            double& resig, double& imsig,
            const double wpscl[], const double ampfac[]);

// --- Internal functions (exposed for use by csigz) ---

/// Single-pole energy-dependent self-energy via adaptive Gaussian quadrature.
///
/// @param ck      Complex momentum
/// @param energy  Energy (complex)
/// @param wi      Plasmon pole energy
/// @param gamma   Broadening (set to 0 in current code)
/// @param amp     Amplitude of plasmon pole
/// @param kfermi  Fermi momentum
/// @param efermi  Fermi energy
FeffComplex sigma1(FeffComplex ck, FeffComplex energy,
                   double wi, double gamma, double amp,
                   double kfermi, double efermi);

/// Derivative of self-energy w.r.t. energy.
FeffComplex dsigma(FeffComplex ck, FeffComplex energy,
                   double wi, double gamma, double amp,
                   double kfermi, double efermi);

/// Hartree-Fock exchange (Dirac-Hara).
///
/// @param ck_in   Complex momentum (absolute units)
/// @param efermi  Fermi energy
/// @param kfermi  Fermi momentum
FeffComplex hfexc(FeffComplex ck_in, double efermi, double kfermi);

// --- Integration kernels ---

/// Plasmon dispersion relation: fq(q) = sqrt((wp - i*gamma)^2 + 4/3*q^2 + q^4)
FeffComplex fq_dispersion(FeffComplex q, const double dppar[10]);

/// Integration kernel r1 (second integral in eq. 13 of H.L.)
FeffComplex kernel_r1(FeffComplex q, const double dppar[10], const FeffComplex cpar[10]);

/// Integration kernel r2 (first integral in eq. 13 of H.L.)
FeffComplex kernel_r2(FeffComplex q, const double dppar[10], const FeffComplex cpar[10]);

/// Integration kernel r3 (third integral, k < kF)
FeffComplex kernel_r3(FeffComplex q, const double dppar[10], const FeffComplex cpar[10]);

/// Derivative kernel dr1
FeffComplex kernel_dr1(FeffComplex q, const double dppar[10], const FeffComplex cpar[10]);

/// Derivative kernel dr2
FeffComplex kernel_dr2(FeffComplex q, const double dppar[10], const FeffComplex cpar[10]);

/// Derivative kernel dr3
FeffComplex kernel_dr3(FeffComplex q, const double dppar[10], const FeffComplex cpar[10]);

/// Type for integration kernel functions
using KernelFn = FeffComplex (*)(FeffComplex q, const double dppar[10],
                                  const FeffComplex cpar[10]);

/// Adaptive Gaussian quadrature for complex integrands.
///
/// Based on Steve White's rewrite of Mike Teter's integration routine,
/// modified by J. Rehr for complex integration.
///
/// @param fn      Function to integrate
/// @param dppar   Double precision parameters passed to fn
/// @param cpar    Complex parameters passed to fn
/// @param xmin    Lower limit of integration
/// @param xmax    Upper limit of integration
/// @param abr     Absolute tolerable error
/// @param rlr     Relative tolerable error
/// @param nsing   Number of singularities
/// @param xsing   Locations of singularities
/// @param error   [out] Error estimate
/// @param numcal  [out] Number of function evaluations
/// @param maxns   [out] Maximum number of regions used
FeffComplex cgratr(KernelFn fn, const double dppar[10], const FeffComplex cpar[10],
                   FeffComplex xmin, FeffComplex xmax,
                   double abr, double rlr,
                   int nsing, FeffComplex xsing[20],
                   double& error, int& numcal, int& maxns);

} // namespace feff::exch
