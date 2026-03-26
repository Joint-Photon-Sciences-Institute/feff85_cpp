#pragma once

// Recursion method (non-symmetric Lanczos) solver for FMS.
// Converted from: ggrm.f
// Also contains helper functions cdot, vecvec, matvec used by all iterative solvers.

#include <vector>
#include <Eigen/Dense>

#include <feff/dimensions.hpp>
#include "fms_types.hpp"

namespace feff::fms {

/// Solve (1 - T*G0) * x = b using the recursion method (Lanczos).
void ggrm(int nsp, const int* i0, int ipi, int ipf, const int* lipotx,
          Eigen::MatrixXcf& g0,
          const Eigen::VectorXcf& tmatrx_diag,
          const Eigen::VectorXcf& tmatrx_off,
          Eigen::MatrixXcf& g0t,
          std::vector<Eigen::MatrixXcf>& gg,
          float toler1, float toler2, const bool* lcalc, int& msord,
          const FMSData& data);

} // namespace feff::fms
