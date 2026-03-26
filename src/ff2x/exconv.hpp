#pragma once
// Convolution with excitation spectrum.
// Converted from: src/FF2X/exconv.f

namespace feff::ff2x {

/// Convolute xmu(E) with excitation spectrum modeled by:
///   f(e) = s02*delta(e) + theta(e)*exp(-e/ed)*x1 + fp(e)
/// Uses fast incremental method exploiting the exponential kernel.
/// Replaces Fortran subroutine exconv.
///
/// @param omega   Energy grid (nk points)
/// @param nk      Number of energy points
/// @param efermi  Fermi level
/// @param s02     Overlap factor
/// @param erelax  Relaxation energy
/// @param wp      Plasmon frequency
/// @param xmu     Input/output absorption coefficient (modified in place)
void exconv(const double omega[], int nk, double efermi,
            double s02, double erelax, double wp, double xmu[]);

} // namespace feff::ff2x
