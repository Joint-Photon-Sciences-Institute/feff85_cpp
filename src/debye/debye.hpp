#pragma once
// Debye-Waller factor calculations.
// Converted from: src/DEBYE/sigms.f, sigcl.f, sigm3.f, sigte3.f, sigrem.f
//
// Provides quantum (Debye model), classical (Debye model),
// third cumulant (Einstein model), and anharmonic Debye-Waller factors
// for multiple-scattering paths.

#include <feff/types.hpp>

namespace feff::debye {

// ---------------------------------------------------------------------------
// Quantum Debye-Waller factor (correlated Debye model)
// Converted from: sigms.f
// ---------------------------------------------------------------------------

/// Calculate sigma^2 for a multiple-scattering path using the
/// quantum correlated Debye model.
///
/// @param tk       Temperature in K
/// @param thetad   Debye temperature in K
/// @param rs       Average Wigner-Seitz/Norman radius in Bohr
/// @param nlegx    Dimension of rat, iz (max legs)
/// @param nleg     Number of legs in path
/// @param rat      Atomic positions rat[0:nleg][3] in Angstroms
///                 (index 0 = index nleg = central atom)
/// @param iz       Atomic numbers iz[0:nleg]
/// @param sig2     Output: Debye-Waller factor in Angstroms^2
void sigms(double tk, double thetad, double rs, int nlegx, int nleg,
           const double rat[][3], const int iz[], double& sig2);

/// Correlation function for quantum Debye model.
/// c(Ri,Rj) = <xi*xj> in the Debye approximation.
void corrfn(double rij, double& cij, double thetad, double tk,
            int iz1, int iz2, double rsavg);

// ---------------------------------------------------------------------------
// Classical Debye-Waller factor (classical Debye model)
// Converted from: sigcl.f
// ---------------------------------------------------------------------------

/// Calculate sigma^2 using the classical Debye model.
/// Same interface as sigms but uses classical coth(x) -> 1/x limit.
void sigcl(double tk, double thetad, double rs, int nlegx, int nleg,
           const double rat[][3], const int iz[], double& sig2);

/// Correlation function for classical Debye model.
void corrfn2(double rij, double& cij, double thetad, double tk,
             int iz1, int iz2, double rsavg);

// ---------------------------------------------------------------------------
// Third cumulant from Einstein model with Morse potential
// Converted from: sigm3.f
// ---------------------------------------------------------------------------

/// Calculate first and third cumulants using correlated Einstein model
/// with a Morse potential (Nguyen Van Hung & J.J.Rehr PRB 56, 43).
///
/// @param sig1     Output: first cumulant
/// @param sig2     Input: second cumulant (sigma^2)
/// @param sig3     Output: third cumulant
/// @param tk       Temperature in K
/// @param alphat   Thermal expansion coefficient (converted to code units internally)
/// @param thetae   Einstein temperature in K
void sigm3(double& sig1, double sig2, double& sig3,
           double tk, double alphat, double thetae);

// ---------------------------------------------------------------------------
// Third cumulant from Debye model (single scattering only)
// Converted from: sigte3.f
// ---------------------------------------------------------------------------

/// Calculate first and third cumulants for single-scattering paths
/// using the Debye model.
///
/// @param iz1      Atomic number of central atom
/// @param iz2      Atomic number of neighbor
/// @param sig2     Second cumulant (sigma^2, code units)
/// @param alphat   Thermal expansion coefficient
/// @param thetad   Debye temperature in K
/// @param reff     Half path length (single precision, from feff.pad)
/// @param sig1     Output: first cumulant
/// @param sig3     Output: third cumulant
void sigte3(int iz1, int iz2, double sig2, double alphat,
            double thetad, double reff, double& sig1, double& sig3);

// ---------------------------------------------------------------------------
// Romberg integration (used internally by corrfn/corrfn2)
// ---------------------------------------------------------------------------

/// Romberg integration of the Debye integrand over [0,1].
/// Used by corrfn (quantum model).
/// @param b     Output: integral value
/// @param eps   Output: estimated relative error
/// @param n     Output: number of refinement steps
void bingrt(double& b, double& eps, int& n);

/// Romberg integration of the classical Debye integrand over [0,1].
/// Used by corrfn2 (classical model).
void bingrt2(double& b, double& eps, int& n);

} // namespace feff::debye
