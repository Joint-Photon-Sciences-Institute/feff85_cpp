#pragma once
// Read feff.inp and populate FeffInput structure.
// Converted from: src/RDINP/rdinp_l.f (~1000 lines)

#include <feff/feff_input.hpp>
#include <string>

namespace feff::rdinp {

/// Read feff.inp file and populate all fields of FeffInput.
/// Returns nabs (number of absorbers for configuration average).
/// Replaces Fortran program/subroutine rdinp.
int rdinp(FeffInput& inp, const std::string& input_file = "feff.inp");

} // namespace feff::rdinp
