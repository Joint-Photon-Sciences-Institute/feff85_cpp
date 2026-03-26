#pragma once
// FeffPipeline: orchestrates the complete FEFF calculation.
// Replaces the Fortran main program feff.f and the individual ffmodN drivers.
//
// Pipeline stages:
//   1. RDINP  — Parse feff.inp, write JSON config files
//   2. POT    — Calculate atomic potentials (SCF)
//   3. XSPH   — Calculate phase shifts and cross-sections
//   4. PATH   — Find scattering paths
//   5. GENFMT — Calculate effective scattering amplitudes (F-matrix)
//   6. FF2X   — Generate output spectra (chi.dat, xmu.dat, etc.)

#include <string>

namespace feff::pipeline {

/// Control flags for which pipeline stages to run.
struct PipelineConfig {
    bool run_pot = true;       // Stage 1+2: parse input + potentials
    bool run_xsph = true;     // Stage 3: phase shifts
    bool run_path = true;     // Stage 4: path finding
    bool run_genfmt = true;   // Stage 5: F-matrix
    bool run_ff2x = true;     // Stage 6: output spectra
    bool verbose = true;       // Enable progress logging
    std::string input_file = "feff.inp";
};

/// Run the complete FEFF pipeline.
/// Returns 0 on success, non-zero on error.
int run_feff(const PipelineConfig& config = PipelineConfig{});

/// Run individual pipeline stages (for modular use).
int run_rdinp(const std::string& input_file = "feff.inp");
int run_pot();
int run_xsph();
int run_path();
int run_genfmt();
int run_ff2x();

} // namespace feff::pipeline
