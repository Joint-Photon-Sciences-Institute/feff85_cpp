#pragma once
// Angular coefficient and radial integral functions for Coulomb/exchange terms.
// Forward declarations for: akeato.f (akeato+bkeato), fdrirk.f
//
// These functions are not yet fully converted; this header provides
// the interface declarations needed by lagdat and other callers.

#include "atom_types.hpp"

namespace feff::atom {

/// Angular coefficient for the direct Coulomb integral F_k(i,j).
/// @param i  First orbital index (1-based)
/// @param j  Second orbital index (1-based)
/// @param k  Multipole order
/// @param angular  Angular coefficients afgk
/// @return  The angular coefficient
double akeato(int i, int j, int k, const AngularCoefficients& angular);

/// Angular coefficient for the exchange Coulomb integral G_k(i,j).
/// @param i  First orbital index (1-based)
/// @param j  Second orbital index (1-based)
/// @param k  Multipole order
/// @param angular  Angular coefficients afgk
/// @return  The angular coefficient (0 if i==j)
double bkeato(int i, int j, int k, const AngularCoefficients& angular);

/// Calculate radial integral R_k(i,j,l,m).
/// R_k = integral of f(r) * u_k(r,s) * g(s)
/// where u_k(r,s) = r_inf^k / r_sup^(k+1)
///
/// @param i,j,l,m  Orbital indices (1-based)
/// @param k        Multipole order
/// @param state    Aggregate atom state
/// @return  Value of the radial integral
double fdrirk(int i, int j, int l, int m, int k, AtomState& state);

} // namespace feff::atom
