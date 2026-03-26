#pragma once
// EXAFS chi(k) calculation using multiple-scattering paths expansion.
// Converted from: src/FF2X/ff2chi.f

#include "ff2x_types.hpp"

namespace feff::ff2x {

/// EXAFS chi(k) calculation.
/// Adds contributions from each path with Debye-Waller factors.
/// Writes chi.dat and xmu.dat output files.
/// Replaces Fortran subroutine ff2chi.
void ff2chi(const FF2xParams& p, int iabs);

} // namespace feff::ff2x
