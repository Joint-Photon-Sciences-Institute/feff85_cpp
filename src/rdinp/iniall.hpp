#pragma once
// Initialize all input variables to their defaults.
// Converted from: src/RDINP/iniall.f

#include <feff/feff_input.hpp>

namespace feff::rdinp {

/// Initialize all fields of FeffInput to their default values.
/// Replaces Fortran subroutine iniall.
void iniall(FeffInput& inp);

} // namespace feff::rdinp
