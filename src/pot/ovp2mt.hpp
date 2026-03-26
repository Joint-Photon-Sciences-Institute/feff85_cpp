#pragma once
// Map overlapped potential/density to muffin-tin spheres.
// Converted from: src/POT/ovp2mt.f
//
// Uses the LU-decomposed overlap matrix from movrlp to decompose
// the overlapped potential (or density) at the MT boundary into
// single-site contributions.

#include <feff/dimensions.hpp>
#include <complex>

namespace feff::pot {

/// Number of overlap grid points (parameter novp in Fortran)
inline constexpr int novp = 40;

/// Dimension of the overlap matrix system
inline constexpr int istatx_dim = novp * (nphx + 1) + 1;

/// Map overlapped potential/density to MT spheres.
///
/// @param nph      Number of unique potentials
/// @param vtot     Potential or density array [251*(nphx+1)], may be overwritten
/// @param lrewr    Rewrite flag:
///                   0 = density calculation (vtot not overwritten)
///                   1 = potential calculation, vint estimated
///                   2 = potential calculation, vint fixed
/// @param qtot     Total electron charge of cluster (density mode only)
/// @param ri       Loucks radial grid [251]
/// @param xnatph   Number of atoms per potential type [nphx+1]
/// @param lnear    Boolean flags per potential [nphx+1]
/// @param inrm     Norman mesh indices [nphx+1]
/// @param imt      MT mesh indices [nphx+1]
/// @param rnrm     Norman radii [nphx+1]
/// @param rmt      MT radii [nphx+1]
/// @param cmovp    LU-decomposed overlap matrix [istatx_dim * istatx_dim]
/// @param ipiv     Pivot indices [istatx_dim]
/// @param vint     Interstitial potential/charge (input/output)
/// @param inters   Interstitial calculation mode
void ovp2mt(int nph, double* vtot, int lrewr, double qtot,
            const double* ri, const double* xnatph, const bool* lnear,
            const int* inrm, const int* imt,
            const double* rnrm, const double* rmt,
            std::complex<float>* cmovp, int* ipiv,
            double& vint, int inters);

} // namespace feff::pot
