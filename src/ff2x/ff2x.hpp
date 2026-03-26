#pragma once
// Main FF2X module interface -- EXAFS chi(k), XANES mu(E), anomalous scattering.
// Converted from: src/FF2X/ff2chi.f, ff2xmu.f, ff2afs.f, ff2gen.f, ffmod6.f

#include "ff2x_types.hpp"

namespace feff::ff2x {

/// Generate spectrum filenames for a given spectrum index iip.
/// Replaces the Fortran filename generation blocks in ff2chi/ff2xmu/ff2afs.
struct SpectrumFiles {
    std::string chi_file;   // e.g. "chi.dat" or "chi0N.dat"
    std::string xmu_file;   // e.g. "xmu.dat" or "xmu0N.dat"
    std::string pad_file;   // e.g. "feff.pad" or "feff0N.bin"
    std::string list_file;  // e.g. "list.dat" or "list0N.dat"
};
SpectrumFiles make_spectrum_files(int iip);

/// Determine if spectrum index is a cross-term.
bool is_cross_spectrum(int iip);

/// EXAFS chi(k) calculation using MS Paths expansion.
/// Replaces Fortran subroutine ff2chi.
void ff2chi(const FF2xParams& p, int iabs);

/// XANES mu(E) calculation using FMS+Paths method.
/// Replaces Fortran subroutine ff2xmu.
void ff2xmu(const FF2xParams& p, int iabs);

/// Anomalous scattering f' / DANES calculation using FMS+Paths.
/// Replaces Fortran subroutine ff2afs.
void ff2afs(const FF2xParams& p, int iabs);

/// Main entry point for FF2X module (module 6).
/// Reads input, dispatches to ff2chi/ff2xmu/ff2afs based on ispec.
/// Replaces Fortran program ffmod6.
void run_ff2x();

} // namespace feff::ff2x
