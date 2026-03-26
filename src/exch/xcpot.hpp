#pragma once

// Main exchange-correlation potential driver
// Converted from src/EXCH/xcpot.f
//
// Dispatches to different XC models: rhl, edp, vbh, rhlbp, csigma/csigz
// based on the index parameter (ixc = index % 10, ibp = index / 10).
//
// Contains internal sigma() dispatcher subroutine converted as a
// static helper function within xcpot.cpp.

#include <feff/types.hpp>
#include <feff/dimensions.hpp>

namespace feff::exch {

/// Main XC potential driver.
///
/// Calculates self-energy correction and complex potential v(r).
///
/// @param iph     Unique potential index (debug/label only)
/// @param ie      Energy index (debug/label only)
/// @param index   XC model selector:
///                  0 = Hedin-Lundqvist + const
///                  1 = Dirac-Hara + const
///                  2 = ground state + const
///                  3 = Dirac-Hara + HL imag + const
///                  ixc = index % 10, ibp = index / 10
/// @param lreal   Nonzero for real self energy
/// @param ifirst  [in/out] First-entry flag; set to 0 before first call
/// @param jri     Index of first interstitial point in Loucks r grid (0-based)
/// @param em      Current energy grid point (complex, a.u.)
/// @param xmu     Fermi level (a.u.)
/// @param vtot    Total potential (Coulomb + gs XC), size nrptx
/// @param vvalgs  Total Coulomb + gs XC from valence electrons, size nrptx
/// @param densty  Electron density, size nrptx
/// @param dmag    Density magnetization, size nrptx
/// @param denval  Valence electron density, size nrptx
/// @param eref    [out] Complex energy reference
/// @param v       [out] Complex potential including energy-dep XC, size nrptx
/// @param vval    [out] As above, XC from valence only, size nrptx
/// @param ipl     Many-pole self energy control
/// @param wpcorr  Plasmon pole frequencies, size MxPole
/// @param ampfac  Plasmon pole amplitudes, size MxPole
/// @param vxcrmu  [workspace] Real part of XC at Fermi level, size nrptx
/// @param vxcimu  [workspace] Imag part of XC at Fermi level, size nrptx
/// @param gsrel   [workspace] Ratio of gs XC with/without magnetization, size nrptx
/// @param vvxcrm  [workspace] Real part of valence XC at Fermi level, size nrptx
/// @param vvxcim  [workspace] Imag part of valence XC at Fermi level, size nrptx
/// @param rnrm    Norman radius (a.u.)
void xcpot(int iph, int ie, int index, int lreal, int& ifirst, int jri,
           FeffComplex em, double xmu,
           const double vtot[], const double vvalgs[],
           const double densty[], const double dmag[], const double denval[],
           FeffComplex& eref, FeffComplex v[], FeffComplex vval[],
           int ipl, const double wpcorr[], const double ampfac[],
           double vxcrmu[], double vxcimu[], double gsrel[],
           double vvxcrm[], double vvxcim[], double rnrm);

} // namespace feff::exch
