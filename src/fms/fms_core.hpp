#pragma once

// Core FMS (Full Multiple Scattering) routines.
// Converted from: fmspack.f (subroutines fms, getkts, xclmz, xgllm)
// Original authors: Bruce Ravel, Alexei Ankudinov
//
// These functions compute the full multiple scattering matrix within a
// cluster at a single energy point.

#include <complex>
#include <vector>
#include <Eigen/Dense>

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>

#include "fms_types.hpp"

namespace feff::fms {

/// Compute full multiple scattering within a cluster at a single energy.
///
/// @param lfms     FMS mode (0=extended/crystal, 1=molecule)
/// @param nsp      Number of spin indices (1 or 2)
/// @param ispin    Spin flag
/// @param inclus   Number of atoms in cluster
/// @param npot     Number of unique potentials
/// @param ck       Complex momentum at this energy [nsp]
/// @param lipotx   Max angular momentum per potential [nphx+1]
/// @param xphase   Phase shifts [nspx][2*lx+1][nphx+1] — flat, Fortran-order
///                 xphase(isp, l, iph) with l in [-lx,lx]
/// @param ik       Energy grid index (for messages)
/// @param iverb    Verbosity (0=silent, 1=print grid info)
/// @param minv     Inversion method: 0=LU, 1=BiCGStab, 2=RM, 3=GM, 4=TFQMR
/// @param rdirec   Direct interaction cutoff distance (Bohr)
/// @param toler1   Convergence tolerance (iterative solvers)
/// @param toler2   Sparsity threshold (iterative solvers)
/// @param lcalc    Which l-channels to compute [lx+1]
/// @param[out] gg  FMS Green's function submatrices
///                 gg(is1, is2, ip) for each unique potential ip
///                 Flat Eigen matrix: gg[ip] is nsp*(lx+1)^2 x nsp*(lx+1)^2
/// @param data     FMS shared data (cluster, rotation, lnlm, dw, basis, cg)
void fms(int lfms, int nsp, int ispin, int inclus, int npot,
         const Complexf* ck, const int* lipotx, const Complexf* xphase,
         int ik, int iverb, int minv, float rdirec,
         float toler1, float toler2, const bool* lcalc,
         std::vector<Eigen::MatrixXcf>& gg,
         FMSData& data);

/// Construct state kets |iat, l, m, isp> for the cluster.
/// Output goes to data.basis.
///
/// @param nsp      Number of spins
/// @param nat      Number of atoms in cluster
/// @param lipotx   Max l per potential [nphx+1]
/// @param[out] i0  Index offset for each potential's representative atom [nphx+1]
/// @param data     FMS shared data (cluster.iphx read, basis written)
void getkts(int nsp, int nat, const int* lipotx, int* i0, FMSData& data);

/// Compute Hankel-like polynomials c_lm(z) by recursion (Rehr-Albers eq.4).
///
/// @param lmaxp1   Largest angular momentum + 1
/// @param mmaxp1   Largest m + 1
/// @param rho      ck * r (complex distance)
/// @param[out] clm Output polynomials [(lx+2) x (2*lx+3)]
void xclmz(int lmaxp1, int mmaxp1, Complexf rho,
            Complexf* clm);

/// Compute G_ll'^|mu|(z) propagator element (Rehr-Albers eq.11,12).
///
/// @param mu       |magnetic state| in sum (mu >= 0)
/// @param ist1     State index 1 (0-based)
/// @param ist2     State index 2 (0-based)
/// @param xclm     c_lm array [(lx+1) x (lx+1) x nclusx x nclusx]
/// @param data     FMS data (basis states, cluster, lnlm)
/// @return         g_ll'^|mu|(z) for this state pair
Complexf xgllm(int mu, int ist1, int ist2,
               const Complexf* xclm, const FMSData& data);

} // namespace feff::fms
