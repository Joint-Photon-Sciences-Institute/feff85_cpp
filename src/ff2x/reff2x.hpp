#pragma once
// Input reading for FF2X module.
// Converted from: src/FF2X/reff2x.f
// Reads ff2x.json and global.json to populate FF2xParams.

#include "ff2x_types.hpp"

namespace feff::ff2x {

/// Read ff2x.json and global.json to populate all FF2X parameters.
/// Replaces Fortran subroutine reff2x.
FF2xParams read_ff2x_input();

} // namespace feff::ff2x
