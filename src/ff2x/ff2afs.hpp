#pragma once
// Anomalous scattering f' / DANES calculation.
// Converted from: src/FF2X/ff2afs.f

#include "ff2x_types.hpp"

namespace feff::ff2x {

/// Anomalous scattering amplitude calculation for a given edge.
/// Uses FMS+Paths method with f' corrections.
/// Writes xmu.dat output.
/// Replaces Fortran subroutine ff2afs.
void ff2afs(const FF2xParams& p, int iabs);

} // namespace feff::ff2x
