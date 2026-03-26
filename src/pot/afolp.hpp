#pragma once
// Automatic overlap fractions.
// Converted from src/POT/afolp.f
//
// Finds folp(iph) automatically and recalculates interstitial parameters.
// Written by ALA 11/97.

#include <feff/dimensions.hpp>

namespace feff::pot {

/// Automatic overlap fraction calculation and parameter update.
///
/// @param verbse  If true, print diagnostic output
/// @param nph     Number of unique potentials
/// @param nat     Number of atoms
/// @param iphat   Potential type per atom [natx]
/// @param rat     Atomic coordinates [3][natx] (column-major)
/// @param iatph   Representative atom per potential [nphx+1]
/// @param xnatph  Atoms per potential type [nphx+1]
/// @param novr    Overlap neighbor counts [nphx+1]
/// @param iphovr  Overlap neighbor pot types [novrx][nphx+1] (flat)
/// @param nnovr   Overlap neighbor multiplicities [novrx][nphx+1] (flat)
/// @param rovr    Overlap neighbor distances [novrx][nphx+1] (flat)
/// @param folp    Overlap fractions [nphx+1] (modified)
/// @param folpx   Max overlap fractions [nphx+1]
/// @param iafolp  Auto-FOLP flag
/// @param edens   Electron density [251][nphx+1] (flat)
/// @param edenvl  Valence density [251][nphx+1] (flat)
/// @param dmag    Spin density [251][nphx+2] (flat)
/// @param vclap   Coulomb potential [251][nphx+1] (flat)
/// @param vtot    Total potential [251][nphx+1] (flat)
/// @param vvalgs  Valence potential [251][nphx+1] (flat)
/// @param imt     MT grid indices [nphx+1]
/// @param inrm    Norman grid indices [nphx+1]
/// @param rmt     MT radii [nphx+1] (modified)
/// @param rnrm    Norman radii [nphx+1]
/// @param ixc     Exchange-correlation model
/// @param rhoint  Interstitial density (output)
/// @param vint    Interstitial potential (output)
/// @param rs      Density parameter (output)
/// @param xf      Fermi momentum (output)
/// @param xmu     Fermi level (input)
/// @param xmunew  New Fermi level (output)
/// @param rnrmav  Average Norman radius (output)
/// @param qtotel  Total charge (output)
/// @param inters  Interstitial model flag
/// @param totvol  Total volume
void afolp(bool verbse, int nph, int nat, const int* iphat, const double* rat,
           const int* iatph, const double* xnatph,
           const int* novr, const int* iphovr, const int* nnovr, const double* rovr,
           double* folp, double* folpx, int iafolp,
           double* edens, double* edenvl,
           double* dmag, double* vclap, double* vtot, double* vvalgs,
           int* imt, int* inrm, double* rmt, double* rnrm,
           int ixc, double& rhoint, double& vint, double& rs, double& xf,
           double xmu, double& xmunew,
           double& rnrmav, double& qtotel, int& inters, double totvol);

} // namespace feff::pot
