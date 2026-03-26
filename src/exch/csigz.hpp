#pragma once

// Non-self-consistent self-energy with renormalization factor Z
// Converted from src/EXCH/csigz.f
//
// Calculates Sigma(k(Energy), Energy) based on an electron gas model
// of epsilon^{-1}, with renormalization Z = 1/(1 - dSigma/dE).
// Uses sigma1 and hfexc from csigma module.

#include <feff/types.hpp>

namespace feff::exch {

/// Non-self-consistent self-energy with renormalization Z.
///
/// Solves: k0^2 = 2*(Energy-Mu) + SigmaF
///         Sigma0 = Sigma(k0, Energy)
///         dSgdE = derivative w.r.t. E
///         k1^2 = k0^2 - 2*(Sigma0-SigmaF)/(1-dSgdE)
///
/// @param energy  Energy at which to evaluate Sigma (complex, a.u.)
/// @param mu      Fermi energy from FEFF (a.u.)
/// @param rs      Density parameter
/// @param resig   [out] Re[Sigma]
/// @param imsig   [out] Im[Sigma]
/// @param ztot    [out] Renormalization factor Z = 1/(1-dSgdE)
/// @param wpscl   Plasmon pole frequency scale factors (MxPole array)
/// @param ampfac  Plasmon pole amplitude factors (MxPole array)
void csigz(FeffComplex energy, double mu, double rs,
           double& resig, double& imsig, FeffComplex& ztot,
           const double wpscl[], const double ampfac[]);

} // namespace feff::exch
