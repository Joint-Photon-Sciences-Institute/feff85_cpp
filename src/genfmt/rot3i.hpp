#pragma once

// Rotation matrix calculation.
// Converted from GENFMT/rot3i.f
//
// Calculates the beta-dependence of rotation matrix elements using
// recursion of formula (4.4.1) in Edmonds.
// First written by J. Mustre (1986), modified by J. Rehr (1987),
// version for genfmt by S. Zabinsky (1991).

#include "genfmt_data.hpp"

namespace feff::genfmt {

/// Calculate rotation matrix elements for leg ileg.
///
/// Input:
///   lxp1    - lmax + 1
///   mxp1    - mmax + 1
///   ileg    - index into beta[] for angle access (1-based from rdpath)
///   beta    - array of beta angles (beta[ileg] is used)
///   istore  - 0-based storage index for rm.dri[...][...][...][istore]
///
/// Output:
///   rm.dri[...][...][...][istore] is set
void rot3i(int lxp1, int mxp1, int ileg, const double beta[],
           RotationMatrixData& rm, int istore);

} // namespace feff::genfmt
