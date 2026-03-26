#pragma once
// Read xsect.bin (cross-section data) into XsectData.
// Converted from: src/FF2X/ff2gen.f (rdxbin subroutine)

#include "ff2x_types.hpp"

namespace feff::ff2x {

/// Read xsect.bin via read_xsect_json and populate XsectData.
/// s02 and mbconv from the caller may override internal values.
/// Replaces Fortran subroutine rdxbin.
XsectData read_xsect_bin(double s02_in, int mbconv);

} // namespace feff::ff2x
