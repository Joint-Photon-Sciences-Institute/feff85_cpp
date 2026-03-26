#pragma once

// Simplified geometric preparation without Debye-Waller factors.
// Converted from: yprep.f
// Used by SCF and LDOS where DW factors should not enter.

#include "fms_types.hpp"

namespace feff::fms {

/// Simplified xprep without Debye-Waller factors.
/// Used in SCF loop and LDOS calculations.
///
/// @param iph0     Potential index of the central atom (0 for absorber)
/// @param nat      Number of atoms in the extended cluster
/// @param[out] inclus  Number of atoms in the FMS cluster
/// @param iphat    Potential index for each atom in extended cluster [nat]
/// @param rmax     FMS cluster radius (Bohr)
/// @param rat      Coordinates (Bohr) [nat][3], row-major
/// @param data     FMS shared state (cluster, rotation, lnlm updated; dw zeroed)
void yprep(int iph0, int nat, int& inclus,
           const int* iphat, float rmax, const float* rat,
           FMSData& data);

} // namespace feff::fms
