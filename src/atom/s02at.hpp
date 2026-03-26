#pragma once
// S02 overlap integral calculation.
// Converted from: src/ATOM/s02at.f

#include "atom_types.hpp"

namespace feff::atom {

/// Calculate S02 shakeup amplitude from orbital overlap integrals.
/// Uses determinant formulation for multi-orbital overlap.
///   ihole: core-hole orbital index
///   norb: number of orbitals
///   nk[30]: kappa quantum numbers
///   xnel[30]: occupations
///   ovpint[30][30]: overlap integrals (input)
///   dval: determinant value (output)
/// Replaces Fortran s02at().
void s02at(int ihole, int norb, const int nk[30], const double xnel[30],
           const double ovpint[][30], double& dval);

/// Helper: construct modified overlap matrix with row/column eliminated.
/// Replaces Fortran elimin().
void elimin(double* d1, int n, double* d2);

} // namespace feff::atom
