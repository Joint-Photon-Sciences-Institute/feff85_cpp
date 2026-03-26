#pragma once
// Low-level numerical utilities for the ATOM module.
// Converted from: aprdev.f, dentfa.f, cofcon.f, potslw.f, messer.f, bkmrdf.f

#include "atom_types.hpp"

namespace feff::atom {

/// Coefficient for term of power (l-1) in the product of two polynomials
/// whose coefficients are in arrays a and b.
/// Fortran 1-based index l; C++ caller passes l in [1..ndor].
/// a and b are 0-based arrays of length max_dev.
double aprdev(const double a[], const double b[], int l);

/// Thomas-Fermi model potential for electrons.
/// dr = distance from nucleus, dz = nuclear charge,
/// ch = ionicity = (number of electrons) - dz - 1.
double dentfa(double dr, double dz, double ch);

/// Convergence acceleration for the iterative SCF process.
/// b is the mixing fraction (kept in [0.1, 0.9]), a = 1 - b.
/// p = current error, q = previous error (updated to p on exit).
void cofcon(double& a, double& b, double& p, double& q);

/// Coulomb potential from electron density using 4-point integration.
/// dv = output potential (length np), d = input density,
/// dr = radial mesh, dpas = exponential step, np = number of points.
void potslw(double dv[], const double d[], const double dr[], double dpas, int np);

/// Print error message from ErrorState to the log.
void messer(ErrorState& err);

/// Angular coefficients for the Breit interaction term.
/// i, j = orbital indices (0-based), k = multipole order.
/// Results stored in breit.cmag[] and breit.cret[].
void bkmrdf(int i, int j, int k, OrbitalConfig& config, BreitCoefficients& breit);

} // namespace feff::atom
