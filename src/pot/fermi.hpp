#pragma once
// Fermi level from interstitial density.
// Converted from src/POT/fermi.f
//
// Calculates mu = V_coulomb(int) + V_xc(int) + kf(int)^2 / 2
// Reference: Lee and Beni, Phys. Rev. B15, 2862 (1977), eq 2.13

namespace feff::pot {

/// Calculate Fermi level from interstitial density.
///
/// @param rhoint  Interstitial density (4*pi*rho)
/// @param vint    Interstitial potential (Coulomb + ground state XC) [Hartrees]
/// @param xmu     Output: Fermi level [Hartrees]
/// @param rs      Output: density parameter
/// @param xf      Output: interstitial Fermi momentum
void fermi(double rhoint, double vint, double& xmu, double& rs, double& xf);

} // namespace feff::pot
