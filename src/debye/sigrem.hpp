#pragma once
// Equation-of-Motion (EM) and Recursion Method (RM) Debye-Waller factors.
// Converted from: src/DEBYE/sigrem.f
// These are advanced DW methods that require spring.inp and feff.inp.
//
// References:
//   EM: Phys. Rev. B, 59, p.948, 1999 (A. Poiarkova)
//   RM: J. Synchrotron Rad., 1999 (A. Poiarkova)

#include <feff/dimensions.hpp>

namespace feff::debye {

// Constants matching Fortran dwpar.h
inline constexpr int natxdw = 200;
inline constexpr int nlegx1 = 9;   // must match legtot
inline constexpr int nphx1  = 7;

/// Equation-of-Motion DW factor calculation.
/// Reads feff.inp and spring.inp on first call, builds dynamical matrix,
/// then solves equations of motion for each path's projected VDOS.
///
/// @param sig2mx   Maximum DW factor seen (updated)
/// @param sig2x    Matrix of max DW by potential pair (updated)
/// @param iem      Unit number for s2_em.dat output (ignored in C++)
/// @param tk       Temperature in K
/// @param ipath    Path index (0 for initialization only)
/// @param nleg     Number of legs
/// @param rat      Path atom positions rat[0:nleg][3] in Angstroms
/// @param sig2     Output: DW factor in Angstroms^2
void sigem(double& sig2mx, double sig2x[],
           int iem, double tk, int ipath, int nleg,
           const double rat[][3], double& sig2);

/// Recursion Method DW factor calculation.
/// Similar interface to sigem but uses recursion method.
void sigrm(double& sig2mx, double sig2x[],
           int ir1, int ir2, double tk, int ipath, int nleg,
           const double rat[][3], double& sig2);

} // namespace feff::debye
