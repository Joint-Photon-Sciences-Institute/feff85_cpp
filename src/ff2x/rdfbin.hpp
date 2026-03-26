#pragma once
// Read feff.pad (Packed ASCII Data format) into FeffPadData.
// Converted from: src/FF2X/rdfbin.f

#include "ff2x_types.hpp"
#include <string>

namespace feff::ff2x {

/// Read feff.pad (or named variant like feff0N.bin).
/// Replaces Fortran subroutine rdfbin.
FeffPadData read_feff_pad(const std::string& filename = "feff.pad");

} // namespace feff::ff2x
