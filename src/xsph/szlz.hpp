#pragma once
// Spin/orbital moment expectation values.
// Converted from src/XSPH/szlz.f
//
// Calculates S_z, L_z, T_z or N_l, N_j-, N_j+ depending on ispin.

#include <feff/types.hpp>
#include <feff/dimensions.hpp>

namespace feff::xsph {

/// Calculate spin/orbital moments and density of states.
///
/// @param verbose  Enable log messages
/// @param ispin   Spin flag: 0=occupation, +-1=Sz/Lz/Tz, +-2=Sz/Nl/Nj
/// @param ecv     Core-valence separation energy
/// @param nph     Number of unique potentials
/// @param nat     Number of atoms
/// @param rgrd    Phase grid step
/// @param nohole  No-hole flag
/// @param rfms2   FMS radius
/// @param lfms2   FMS angular momentum cutoff
/// @param lmaxph  Max l per potential [nphx+1]
/// @param edens   Electron density [251][nphx+1]
/// @param edenvl  Valence density [251][nphx+1]
/// @param dmag    Magnetization density [251][nphx+2]
/// @param vtot    Total potential [251][nphx+1]
/// @param vvalgs  Valence potential [251][nphx+1]
/// @param rmt     Muffin-tin radii [nphx+1]
/// @param rnrm    Norman radii [nphx+1]
/// @param ixc     Exchange-correlation model
/// @param rhoint  Interstitial density
/// @param vint    Interstitial potential
/// @param xmu     Fermi level (Hartrees)
/// @param jumprm  Jump flag for potential
/// @param xnval   Valence occupations [30][nphx+1]
/// @param iorb    Orbital indices [-4:3][nphx+1]
/// @param x0      Grid parameter
/// @param dx      Original grid step
/// @param xion    Ionicities [nphx+1]
/// @param iunf    Unfolding flag
/// @param iz      Atomic numbers [nphx+1]
/// @param adgc    Development coefficients [10][30][nphx+2]
/// @param adpc    Development coefficients [10][30][nphx+2]
/// @param dgc     Dirac components [251][30][nphx+2]
/// @param dpc     Dirac components [251][30][nphx+2]
/// @param ihole   Core-hole index
/// @param rat     Atomic coordinates [3][natx]
/// @param iphat   Potential assignments [natx]
/// @param corr    Output: correction factor for spin-orbit
void szlz(bool verbose, int ispin, double ecv, int nph, int nat,
          double rgrd, int nohole, float rfms2, int lfms2,
          const int lmaxph[], double* edens, double* edenvl,
          double* dmag, double* vtot, double* vvalgs,
          const double rmt[], const double rnrm[],
          int ixc, double rhoint, double vint, double xmu, int jumprm,
          const double* xnval, const int* iorb,
          double x0, double dx, const double xion[], int iunf, const int iz[],
          const double* adgc, const double* adpc,
          const double* dgc, const double* dpc,
          int ihole, const double* rat, const int iphat[], double& corr);

} // namespace feff::xsph
