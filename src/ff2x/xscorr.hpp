#pragma once
// Cross-section correction via contour integration.
// Converted from: src/FF2X/xscorr.f
// Convolutes xmu(E) with Lorentzian using complex energy plane calculations.

#include <feff/types.hpp>

namespace feff::ff2x {

/// Convolute xmu(E) = xsec + xsnorm*chia with Lorentzian using
/// calculations in the complex energy plane.
/// Replaces Fortran subroutine xscorr.
///
/// @param ispec   Spectroscopy type
/// @param emxs    Complex energy grid (ne points)
/// @param ne1     Number of horizontal axis points
/// @param ne      Total number of points
/// @param ik0     Fermi level index
/// @param xsec    Cross-section (ne, modified in place for output)
/// @param xsnorm  Normalized cross-section (ne)
/// @param chia    Chi array (ne)
/// @param vrcorr  Real energy correction
/// @param vicorr  Imaginary energy correction (set to 0 internally)
/// @param cchi    Output: correction array (ne)
void xscorr(int ispec, FeffComplex emxs[], int ne1, int ne, int ik0,
            FeffComplex xsec[], const double xsnorm[], FeffComplex chia[],
            double vrcorr, double vicorr, FeffComplex cchi[]);

/// Lorentzian function used in xscorr and fprime.
///   lorenz(xloss, w, dele) = xloss / pi / (xloss^2 + (i*w - dele)^2)
FeffComplex lorenz(double xloss, double w, double dele);

/// Step function used in xscorr.
///   astep(xloss, dele) = 0.5 + atan(dele/xloss) / pi
double astep(double xloss, double dele);

} // namespace feff::ff2x
