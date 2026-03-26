#pragma once
// Set radial grid indices for MT and Norman radii.
// Converted from src/POT/sidx.f

namespace feff::pot {

/// Find grid indices for muffin-tin and Norman radii, and the last
/// non-zero density point.
///
/// @param rholap  Overlapped density array [npts]
/// @param npts    Size of rholap
/// @param rmt     Muffin-tin radius (may be modified if density zero inside)
/// @param rnrm    Norman radius (may be modified if density zero inside)
/// @param imax    Output: last non-zero density index (1-based)
/// @param imt     Output: grid index for muffin-tin radius (1-based)
/// @param inrm    Output: grid index for Norman radius (1-based)
void sidx(const double* rholap, int npts, double& rmt, double& rnrm,
          int& imax, int& imt, int& inrm);

} // namespace feff::pot
