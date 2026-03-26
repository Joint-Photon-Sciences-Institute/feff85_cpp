#pragma once

// LU decomposition solver for the FMS Green's function matrix.
// Converted from: gglu.f
// Uses Eigen PartialPivLU instead of LAPACK cgetrf/cgetrs.

#include <vector>
#include <Eigen/Dense>

#include <feff/dimensions.hpp>
#include "fms_types.hpp"

namespace feff::fms {

/// Solve (1 - G0*T) * Gc = G0 by LU decomposition.
///
/// @param nsp          Number of spins
/// @param i0           Index offset per potential [nphx+1]
/// @param ipi, ipf     Potential index range to compute
/// @param lipotx       Max l per potential [nphx+1]
/// @param g0           Free propagator matrix [istate x istate]
/// @param tmatrx_diag  Diagonal T-matrix elements [istate]
/// @param tmatrx_off   Off-diagonal T-matrix elements [istate] (spin-flip)
/// @param g0t          Work matrix [istate x istate] (overwritten)
/// @param[out] gg      Output Green's function per potential
/// @param data         FMS data (basis states)
void gglu(int nsp, const int* i0, int ipi, int ipf, const int* lipotx,
          Eigen::MatrixXcf& g0,
          const Eigen::VectorXcf& tmatrx_diag,
          const Eigen::VectorXcf& tmatrx_off,
          Eigen::MatrixXcf& g0t,
          std::vector<Eigen::MatrixXcf>& gg,
          const FMSData& data);

} // namespace feff::fms
