#pragma once

// Full FMS orchestrator — drives the energy loop for XANES/FMS calculations.
// Converted from: fmstot.f
// Original author: Alexei Ankudinov (06.1997)
//
// This reads phase.pad, prepares the cluster geometry, and loops over
// all energy points calling the FMS solver. Results are written to fms.bin.

#include <complex>
#include <vector>
#include <array>

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <feff/types.hpp>

#include "fms_types.hpp"

namespace feff::fms {

/// Run the full FMS calculation over all energy points.
///
/// Reads phase.pad, constructs the cluster, calls fms() for each energy,
/// and writes fms.bin.
///
/// @param rclust   FMS cluster radius (Angstrom, converted to Bohr internally)
/// @param idwopt   Debye-Waller option
/// @param tk       Temperature (K)
/// @param thetad   Debye temperature (K)
/// @param sigma2   Global sigma^2 (Ang^2)
/// @param lmaxph   Max angular momentum per potential [nphx+1]
/// @param nat      Number of atoms
/// @param iphat    Potential type per atom [nat]
/// @param ratdbl   Coordinates (Bohr, double) [nat][3]
/// @param ipol     Polarization type
/// @param ispin    Spin flag
/// @param le2      Multipole flag
/// @param angks    Angle between k-vector and spin-vector
/// @param ptz      Polarization tensor [3][3]
/// @param minv     Matrix inversion method
/// @param rdirec   Direct interaction cutoff (Angstrom)
/// @param toler1   Convergence tolerance
/// @param toler2   Sparsity threshold
void fmstot(float rclust, int idwopt, double tk, double thetad, double sigma2,
            int* lmaxph, int nat, const int* iphat, const double* ratdbl,
            int ipol, int ispin, int le2, double angks,
            const FeffComplex ptz[3][3],
            int minv, float rdirec, float toler1, float toler2);

} // namespace feff::fms
