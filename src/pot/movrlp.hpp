#pragma once
// Overlap matrix construction for muffin-tin potential decomposition.
// Converted from: src/POT/movrlp.f
//
// Constructs the overlap matrix based on geometry of overlapped
// muffin-tin spheres and performs LU decomposition for use by ovp2mt.

#include <feff/dimensions.hpp>
#include "ovp2mt.hpp"  // for novp, istatx_dim
#include <complex>

namespace feff::pot {

/// Construct and LU-decompose the overlap matrix.
///
/// @param nph      Number of unique potentials
/// @param nat      Number of atoms in cluster
/// @param iphat    Potential index per atom [natx]
/// @param rat      Atomic coordinates [3*natx], column-major
/// @param iatph    Representative atom per potential [nphx+1]
/// @param xnatph   Stoichiometry per potential type [nphx+1]
/// @param novr     Number of overlap shells per potential [nphx+1]
/// @param iphovr   Overlap shell potential indices [novrx*(nphx+1)]
/// @param nnovr    Overlap shell atom counts [novrx*(nphx+1)]
/// @param rovr     Overlap shell radii [novrx*(nphx+1)]
/// @param imt      MT mesh indices [nphx+1]
/// @param rmt      MT radii [nphx+1]
/// @param rnrm     Norman radii [nphx+1]
/// @param ri       Radial grid [251] (output, filled with Loucks grid)
/// @param lnear    Per-potential boolean flags [nphx+1] (output)
/// @param cmovp    Overlap matrix [istatx_dim*istatx_dim] (output, LU-decomposed)
/// @param ipiv     Pivot indices [istatx_dim] (output)
/// @param volint   Interstitial volume correction (input/output)
/// @param inters   Interstitial calculation mode
void movrlp(int nph, int nat, const int* iphat, const double* rat,
            const int* iatph, const double* xnatph,
            const int* novr, const int* iphovr, const int* nnovr,
            const double* rovr,
            const int* imt, const double* rmt, const double* rnrm,
            double* ri, bool* lnear,
            std::complex<float>* cmovp, int* ipiv,
            double& volint, int inters);

} // namespace feff::pot
