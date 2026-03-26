#pragma once
// Convert experimental units (Angstrom, eV) to atomic units (bohr, hartree).
// Converted from: src/RDINP/eu2au.f

#include <feff/feff_input.hpp>

namespace feff::rdinp {

/// Convert all relevant fields of FeffInput from experimental units
/// (Angstrom, eV) to atomic units (bohr, hartree).
/// Note: sig2g is intentionally NOT transformed.
/// Note: this operates on the sorted geometry (rat/nat), not on ratx/natt.
///       The caller must pass rat/nat separately since they are local to ffsort.
/// Replaces Fortran subroutine eu2au.
void eu2au(FeffInput& inp, int nat, double rat[][3]);

} // namespace feff::rdinp
