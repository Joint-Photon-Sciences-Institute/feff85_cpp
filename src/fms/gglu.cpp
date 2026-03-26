// LU decomposition solver for FMS Green's function.
// Converted from: gglu.f
// Uses Eigen PartialPivLU instead of LAPACK cgetrf/cgetrs.

#include "gglu.hpp"

#include <stdexcept>
#include <iostream>

namespace feff::fms {

void gglu(int nsp, const int* i0, int ipi, int ipf, const int* lipotx,
          Eigen::MatrixXcf& g0,
          const Eigen::VectorXcf& tmatrx_diag,
          const Eigen::VectorXcf& tmatrx_off,
          Eigen::MatrixXcf& g0t,
          std::vector<Eigen::MatrixXcf>& gg,
          const FMSData& data) {

    const auto& basis = data.basis;
    int istate = basis.istate;

    // Construct g0t = I - G0 * T
    // T is (nearly) diagonal: tmatrx_diag on diagonal, tmatrx_off for spin-flip
    for (int icol = 0; icol < istate; ++icol) {
        for (int irow = 0; irow < istate; ++irow) {
            // Diagonal T contribution
            g0t(irow, icol) = -g0(irow, icol) * tmatrx_diag(icol);

            // Off-diagonal T contribution (spin-flip, nsp==2 only)
            int l1   = basis.l(icol);
            int m1   = basis.m(icol);
            int isp1 = basis.isp(icol);
            int m2   = m1 + isp1;
            if (nsp == 2 && m2 > -l1 + 1 && m2 < l1 + 2) {
                // ist2 = icol + (-1)^(isp1+1)  (Fortran 1-based isp, C++ 0-based)
                int ist2 = icol + ((isp1 == 0) ? 1 : -1);
                if (ist2 >= 0 && ist2 < istate) {
                    g0t(irow, icol) -= g0(irow, ist2) * tmatrx_off(icol);
                }
            }
        }
        g0t(icol, icol) += Complexf(1.0f, 0.0f);
    }

    // LU decomposition using Eigen
    Eigen::PartialPivLU<Eigen::MatrixXcf> lu(g0t.topLeftCorner(istate, istate));

    // For each potential, solve g0t * result = g0_columns and extract submatrix
    for (int ip = ipi; ip <= ipf; ++ip) {
        int ipart = nsp * (lipotx[ip] + 1) * (lipotx[ip] + 1);
        int offset = i0[ip];

        // Extract the relevant columns of g0
        Eigen::MatrixXcf g0_cols = g0.block(0, offset, istate, ipart);

        // Solve: g0t * x = g0_cols  =>  x = g0t^{-1} * g0_cols
        Eigen::MatrixXcf result = lu.solve(g0_cols);

        // Pack into gg: extract rows [offset, offset+ipart)
        for (int is2 = 0; is2 < ipart; ++is2) {
            for (int is1 = 0; is1 < ipart; ++is1) {
                gg[ip](is1, is2) = result(is1 + offset, is2);
            }
        }
    }
}

} // namespace feff::fms
