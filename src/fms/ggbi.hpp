#pragma once

// BiCGStab (Bi-Conjugate Gradient Stabilized) solver for FMS.
// Converted from: ggbi.f
// Reference: Saad, "Iterative Methods for Sparse Matrices", p.220 (1996)

#include <vector>
#include <Eigen/Dense>

#include <feff/dimensions.hpp>
#include "fms_types.hpp"

namespace feff::fms {

/// Solve (1 - T*G0) * x = b using BiCGStab iteration.
///
/// @param nsp          Number of spins
/// @param i0           Index offset per potential [nphx+1]
/// @param ipi, ipf     Potential index range
/// @param lipotx       Max l per potential [nphx+1]
/// @param g0           Free propagator [istate x istate]
/// @param tmatrx_diag  Diagonal T-matrix [istate]
/// @param tmatrx_off   Off-diagonal T-matrix [istate]
/// @param g0t          Work matrix (1 - T*G0) [istate x istate]
/// @param[out] gg      Output Green's function per potential
/// @param toler1       Convergence tolerance
/// @param toler2       Sparsity threshold
/// @param lcalc        Which l-channels to compute [lx+1]
/// @param[out] msord   Total number of matrix-vector products
/// @param data         FMS data
void ggbi(int nsp, const int* i0, int ipi, int ipf, const int* lipotx,
          Eigen::MatrixXcf& g0,
          const Eigen::VectorXcf& tmatrx_diag,
          const Eigen::VectorXcf& tmatrx_off,
          Eigen::MatrixXcf& g0t,
          std::vector<Eigen::MatrixXcf>& gg,
          float toler1, float toler2, const bool* lcalc, int& msord,
          const FMSData& data);

} // namespace feff::fms
