#pragma once
// F-prime (anomalous scattering) calculation.
// Converted from: src/FF2X/fprime.f
// Calculates f' including solid state and lifetime effects
// using algorithm in Ankudinov, Rehr DANES paper.

#include <feff/types.hpp>

namespace feff::ff2x {

/// Calculate f' including solid state and lifetime effects.
/// Output correction is returned via cchi. The rest is input.
///   mu(omega) = xsec + xsnorm*chia + cchi
///
/// @param ei      Fermi energy (code units)
/// @param emxs    Complex energy grid (ne points)
/// @param ne1     Number of horizontal axis points
/// @param ne3     Number of auxiliary axis points
/// @param ne      Total number of energy points
/// @param ik0     Fermi level index
/// @param xsec    Cross-section array (ne, may be modified)
/// @param xsnorm  Normalized cross-section (ne)
/// @param chia    Chi array (ne)
/// @param vrcorr  Real energy correction
/// @param vicorr  Imaginary energy correction
/// @param cchi    Output: correction array (ne)
void fprime(double ei, FeffComplex emxs[], int ne1, int ne3, int ne, int ik0,
            FeffComplex xsec[], double xsnorm[], FeffComplex chia[],
            double vrcorr, double vicorr, FeffComplex cchi[]);

/// Anomalous f' logarithmic function.
/// @param icase  1=simplified, 2=real w, 3=pure imaginary w
FeffComplex funlog(int icase, double xloss, double w, double dele);

/// Integration for f' between points n1 and n2 on vertical axis.
void fpint(const FeffComplex emxs[], const FeffComplex xmu[],
           int n1, int n2, double dele, double xloss, double eps4_val,
           double efermi, FeffComplex& value);

/// Integration for f' between points 1 and n2 on horizontal axis,
/// plus tail to infinity.
void fpintp(const double em[], const FeffComplex xmu[],
            int n2, double dele, double xloss, double efermi,
            FeffComplex& value);

} // namespace feff::ff2x
