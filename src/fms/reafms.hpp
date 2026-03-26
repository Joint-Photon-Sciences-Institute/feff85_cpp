#pragma once

// Read FMS configuration from fms.json.
// Converted from: reafms.f
// Uses nlohmann/json (matching the project's JSON I/O pattern).

#include <complex>
#include <array>

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <feff/types.hpp>

namespace feff::fms {

/// FMS configuration parameters read from fms.json.
struct FmsConfig {
    int mfms    = 0;     // FMS module flag
    int idwopt  = 0;     // Debye-Waller option
    int minv    = 0;     // Matrix inversion method

    float rfms2   = 0.0f;  // FMS radius (Angstrom, converted to Bohr by caller)
    float rdirec  = 0.0f;  // Direct cutoff (Angstrom, converted to Bohr by caller)
    float toler1  = 0.0f;  // Convergence tolerance
    float toler2  = 0.0f;  // Sparsity threshold

    double tk     = 0.0;   // Temperature (K)
    double thetad = 0.0;   // Debye temperature (K)
    double sig2g  = 0.0;   // Global sigma^2 (Ang^2)

    std::array<int, nphx + 1> lmaxph{};  // Max angular momentum per potential
};

/// Read FMS configuration from fms.json, geometry from geom.json,
/// and global parameters from global.json.
///
/// Also reads: nat, iphat, rat, ipol, ispin, le2, angks, ptz
/// (same as Fortran reafms).
///
/// @param[out] config   FMS configuration
/// @param[out] nat      Number of atoms
/// @param[out] iphat    Potential type per atom (caller-allocated, size natx)
/// @param[out] rat      Coordinates (Angstrom, double) [natx][3]
/// @param[out] ipol     Polarization type
/// @param[out] ispin    Spin flag
/// @param[out] le2      Multipole flag
/// @param[out] angks    Angle between k and spin vectors
/// @param[out] ptz      Polarization tensor [3][3]
void reafms(FmsConfig& config,
            int& nat, int* iphat, double* rat,
            int& ipol, int& ispin, int& le2, double& angks,
            std::array<std::array<FeffComplex, 3>, 3>& ptz);

} // namespace feff::fms
