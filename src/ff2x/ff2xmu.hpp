#pragma once
// XANES mu(E) calculation using FMS+Paths method.
// Converted from: src/FF2X/ff2xmu.f

#include "ff2x_types.hpp"

namespace feff::ff2x {

/// XANES mu(E) calculation using FMS+Paths.
/// Uses the FMS Green's function (gtr from fms.bin) plus path contributions.
/// Writes xmu.dat output.
/// Replaces Fortran subroutine ff2xmu.
void ff2xmu(const FF2xParams& p, int iabs);

} // namespace feff::ff2x
