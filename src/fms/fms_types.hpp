#pragma once

// FMS module types — structs replacing Fortran COMMON blocks.
// Converted from: xstruc.h, xparam.h, and common blocks in fmspack.f / xprep.f
//
// All arrays use 0-based indexing.  The Fortran code used 1-based atoms and
// various negative-index dimensions; equivalent C++ layouts are flat or
// use Eigen matrices.

#include <complex>
#include <vector>
#include <array>
#include <Eigen/Dense>

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>

namespace feff::fms {

// Single-precision complex used throughout FMS (Fortran: complex)
using Complexf = std::complex<float>;

// Derived dimension constants (from xparam.h)
inline constexpr int nphasx = nphx;
inline constexpr int natxx  = natx;
inline constexpr int nexx   = nex;
inline constexpr int istatx = nspx * (lx + 1) * (lx + 1) * nclusx;

// ---------------------------------------------------------------------------
// ClusterData — replaces COMMON /xstruc/
//   xphi(nclusx, nclusx) — angles between z-axis and atom pairs
//   xrat(3, nclusx)      — Cartesian coordinates in Bohr
//   iphx(nclusx)         — potential index per atom
// ---------------------------------------------------------------------------
struct ClusterData {
    // xphi[i][j]: angle between z-axis and vector R_i - R_j
    // Stored row-major: xphi[i * nclusx + j]
    std::vector<float> xphi;   // nclusx * nclusx

    // xrat[i][dim]: coordinate dim (0=x,1=y,2=z) of atom i
    // Stored as xrat[i * 3 + dim]
    std::vector<float> xrat;   // nclusx * 3

    // iphx[i]: unique-potential index for atom i
    std::vector<int> iphx;     // nclusx

    ClusterData()
        : xphi(nclusx * nclusx, 0.0f),
          xrat(nclusx * 3, 0.0f),
          iphx(nclusx, 0)
    {}

    // Accessors with bounds notation matching Fortran (0-based atom index)
    float& phi(int i, int j)       { return xphi[i * nclusx + j]; }
    float  phi(int i, int j) const { return xphi[i * nclusx + j]; }

    float& rat(int iat, int dim)       { return xrat[iat * 3 + dim]; }
    float  rat(int iat, int dim) const { return xrat[iat * 3 + dim]; }
};

// ---------------------------------------------------------------------------
// RotationData — replaces COMMON /rotx/ and COMMON /rotsav/
//   drix(-lx:lx, -lx:lx, 0:lx, 0:1, nclusx, nclusx)
// ---------------------------------------------------------------------------
struct RotationData {
    // Flat storage for 6D rotation matrix array.
    // Dimensions: (2*lx+1) x (2*lx+1) x (lx+1) x 2 x nclusx x nclusx
    static constexpr int dim_m  = 2 * lx + 1;
    static constexpr int dim_l  = lx + 1;
    static constexpr int dim_k  = 2;
    static constexpr int stride = dim_m * dim_m * dim_l * dim_k;
    std::vector<Complexf> drix;  // stride * nclusx * nclusx

    // Saved rotation matrices cache (replaces /rotsav/)
    static constexpr int jsavx = 150;
    static constexpr float roteps = 1.0e-12f;
    static constexpr int jbmagk = -9999;

    std::vector<float> drisav;   // dim_m * dim_m * dim_l * jsavx
    std::vector<float> betsav;   // jsavx
    std::vector<int> ldsav;      // jsavx
    std::vector<int> mdsav;      // jsavx
    int jsav = 0;

    RotationData()
        : drix(static_cast<size_t>(stride) * nclusx * nclusx, Complexf(0, 0)),
          drisav(static_cast<size_t>(dim_m) * dim_m * dim_l * jsavx, 0.0f),
          betsav(jsavx, static_cast<float>(jbmagk)),
          ldsav(jsavx, 0),
          mdsav(jsavx, 0)
    {}

    // Access drix[m2+lx][m1+lx][il][k][j][i]
    // Fortran: drix(m2, m1, il, k, j, i) with m2,m1 in [-lx,lx], il in [0,lx], k in [0,1]
    int drix_index(int m2, int m1, int il, int k, int j, int i) const {
        return (m2 + lx)
             + dim_m * ((m1 + lx)
             + dim_m * (il
             + dim_l * (k
             + dim_k * (j
             + nclusx * i))));
    }

    Complexf& at(int m2, int m1, int il, int k, int j, int i) {
        return drix[drix_index(m2, m1, il, k, j, i)];
    }
    Complexf at(int m2, int m1, int il, int k, int j, int i) const {
        return drix[drix_index(m2, m1, il, k, j, i)];
    }

    // Access saved rotation matrices
    int drisav_index(int m2, int m1, int il, int isav) const {
        return (m2 + lx) + dim_m * ((m1 + lx) + dim_m * (il + dim_l * isav));
    }
};

// ---------------------------------------------------------------------------
// LegendreNorm — replaces COMMON /lnlm/
//   xnlm(0:lx, 0:lx) — normalization factors
// ---------------------------------------------------------------------------
struct LegendreNorm {
    std::array<float, (lx + 1) * (lx + 1)> xnlm{};

    float& at(int m, int l)       { return xnlm[m + (lx + 1) * l]; }
    float  at(int m, int l) const { return xnlm[m + (lx + 1) * l]; }
};

// ---------------------------------------------------------------------------
// DebyeWaller — replaces COMMON /xdwf/
//   sigsqr(nclusx, nclusx) — pair-wise mean-square displacements (Ang^2)
// ---------------------------------------------------------------------------
struct DebyeWaller {
    std::vector<float> sigsqr;  // nclusx * nclusx

    DebyeWaller() : sigsqr(nclusx * nclusx, 0.0f) {}

    float& at(int i, int j)       { return sigsqr[i * nclusx + j]; }
    float  at(int i, int j) const { return sigsqr[i * nclusx + j]; }
};

// ---------------------------------------------------------------------------
// BasisStates — replaces COMMON /stkets/
//   lrstat(4, istatx) — state kets |iat, l, m, isp>
//   istate             — number of active states
// ---------------------------------------------------------------------------
struct BasisStates {
    // lrstat[ist][0..3] = {iat, l, m, isp}  (0-based state index)
    std::vector<std::array<int, 4>> lrstat;
    int istate = 0;

    BasisStates() : lrstat(istatx) {}

    // Convenience accessors (0-based state index)
    int iat(int ist) const { return lrstat[ist][0]; }
    int l  (int ist) const { return lrstat[ist][1]; }
    int m  (int ist) const { return lrstat[ist][2]; }
    int isp(int ist) const { return lrstat[ist][3]; }
};

// ---------------------------------------------------------------------------
// ClebschGordon — replaces COMMON /t3j/
//   t3jp(0:lx, -lx:lx, 2) — <LS|J> for j = l+1/2
//   t3jm(0:lx, -lx:lx, 2) — <LS|J> for j = l-1/2
// ---------------------------------------------------------------------------
struct ClebschGordon {
    // Flat storage: [l][m+lx][isp]
    static constexpr int dim_l = lx + 1;
    static constexpr int dim_m = 2 * lx + 1;

    std::array<float, dim_l * dim_m * 2> t3jp{};
    std::array<float, dim_l * dim_m * 2> t3jm{};

    int index(int l, int m, int isp) const {
        return l + dim_l * ((m + lx) + dim_m * isp);
    }

    float  jp(int l, int m, int isp) const { return t3jp[index(l, m, isp)]; }
    float& jp(int l, int m, int isp)       { return t3jp[index(l, m, isp)]; }
    float  jm(int l, int m, int isp) const { return t3jm[index(l, m, isp)]; }
    float& jm(int l, int m, int isp)       { return t3jm[index(l, m, isp)]; }
};

// ---------------------------------------------------------------------------
// FMSData — aggregation of all shared state for FMS calculation
// ---------------------------------------------------------------------------
struct FMSData {
    ClusterData   cluster;
    RotationData  rotation;
    LegendreNorm  lnlm;
    DebyeWaller   dw;
    BasisStates   basis;
    ClebschGordon cg;
};

} // namespace feff::fms
