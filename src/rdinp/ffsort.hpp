#pragma once
// Sort atoms, build geometry, and create polarization tensor.
// Converted from: src/RDINP/ffsort.f

#include <feff/feff_input.hpp>

namespace feff::rdinp {

/// Find the iabs-th absorber atom in the full atom list, build a sorted
/// geometry of atoms within rclabs of that absorber, and optionally
/// construct the polarization tensor via mkptz.
///
/// @param iabs   Which absorber atom to use (1-based index among iphabs-type atoms)
/// @param doptz  If true, call mkptz to build polarization tensor
/// @param inp    Master input structure (modified: nph, ptz, angks, le2, evec, etc.)
///
/// Replaces Fortran subroutine ffsort(iabs, doptz).
void ffsort(int iabs, bool doptz, FeffInput& inp);

} // namespace feff::rdinp
