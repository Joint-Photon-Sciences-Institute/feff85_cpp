// Core FMS (Full Multiple Scattering) routines.
// Converted from: fmspack.f
// Original authors: Bruce Ravel, Alexei Ankudinov (feb 2000)
//
// Key changes from Fortran:
//   - COMMON blocks replaced by FMSData struct
//   - 1-based indexing converted to 0-based
//   - LAPACK (cgetrf/cgetrs) replaced by Eigen
//   - Large work arrays (g0, g0t) allocated dynamically

#include "fms_core.hpp"
#include "gglu.hpp"
#include "ggbi.hpp"
#include "ggrm.hpp"
#include "gggm.hpp"
#include "ggtf.hpp"

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <iostream>

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>

namespace feff::fms {

// ============================================================================
// getkts — construct state kets |iat, l, m, isp>
// ============================================================================
void getkts(int nsp, int nat, const int* lipotx, int* i0, FMSData& data) {
    auto& basis = data.basis;
    const auto& iphx = data.cluster.iphx;

    basis.istate = 0;

    for (int iat = 0; iat < nat; ++iat) {
        int ip = iphx[iat];
        // i0[ip]: index of the first state for this potential's representative atom
        if (i0[ip] < 0) i0[ip] = basis.istate;
        int lim = std::min(lx, lipotx[ip]);
        for (int l = 0; l <= lim; ++l) {
            for (int m = -l; m <= l; ++m) {
                for (int isp = 0; isp < nsp; ++isp) {
                    if (basis.istate >= istatx) {
                        throw std::runtime_error(
                            "getkts: Exceeded maximum number of LR states (istatx=" +
                            std::to_string(istatx) + ")");
                    }
                    basis.lrstat[basis.istate] = {iat, l, m, isp};
                    ++basis.istate;
                }
            }
        }
    }
}

// ============================================================================
// xclmz — Hankel-like polynomials c_lm(z) by recursion
// Rehr-Albers eq.4.  c(il,im) = c_l^(m) z^m/m!
// Array layout: clm[(lx+2)][(2*lx+3)], stored row-major clm[il * (2*lx+3) + im]
// Indices start at 0 (shifted from Fortran's 1-based).
// ============================================================================
void xclmz(int lmaxp1, int mmaxp1, Complexf rho, Complexf* clm) {
    constexpr int ltotb = lx + 1;
    constexpr int mntot = 2 * ltotb;  // mtotb + ntotb
    constexpr int cols = mntot + 1;
    constexpr int rows = ltotb + 1;

    Complexf coni_f(0.0f, 1.0f);

    // Initialize to zero
    for (int i = 0; i < rows * cols; ++i) clm[i] = Complexf(0.0f, 0.0f);

    auto C = [&](int il, int im) -> Complexf& {
        return clm[il * cols + im];
    };

    Complexf z = (-coni_f) / rho;
    Complexf cmm(1.0f, 0.0f);

    C(0, 0) = Complexf(1.0f, 0.0f);
    C(1, 0) = C(0, 0) - z;

    int lmax = lmaxp1 - 1;
    for (int il = 2; il <= lmax; ++il) {
        C(il, 0) = C(il - 2, 0) - z * static_cast<float>(2 * il - 1) * C(il - 1, 0);
    }

    int mmxp1 = std::min(lmaxp1, mmaxp1);
    for (int im = 1; im < mmxp1; ++im) {
        cmm = (-cmm) * static_cast<float>(2 * im - 1) * z;
        C(im, im) = cmm;
        C(im + 1, im) = cmm * static_cast<float>(2 * im + 1) *
                         (1.0f - static_cast<float>(im + 1) * z);
        for (int il = im + 1; il <= lmax; ++il) {
            C(il, im) = C(il - 2, im) - static_cast<float>(2 * il - 1) * z *
                         (C(il - 1, im) + C(il - 1, im - 1));
        }
    }
}

// ============================================================================
// xgllm — G_ll'^|mu|(z) from Rehr-Albers eq.11,12
// ============================================================================
Complexf xgllm(int mu, int ist1, int ist2,
               const Complexf* xclm, const FMSData& data) {
    // xclm layout: [m][l][j][i] each in appropriate ranges
    // Flattened as xclm[m + (lx+1) * (l + (lx+1) * (j + nclusx * i))]
    auto xclm_at = [&](int m, int l, int j, int i) -> Complexf {
        return xclm[m + (lx + 1) * (l + (lx + 1) * (j + nclusx * i))];
    };

    const auto& basis = data.basis;
    const auto& lnlm  = data.lnlm;

    int iat1 = basis.iat(ist1);
    int l1   = basis.l(ist1);
    int iat2 = basis.iat(ist2);
    int l2   = basis.l(ist2);
    int numax = std::min(l1, l2 - mu);

    Complexf sum(0.0f, 0.0f);
    for (int nu = 0; nu <= numax; ++nu) {
        int mn = mu + nu;

        Complexf gamtl = static_cast<float>(2 * l1 + 1) *
                         xclm_at(nu, l1, iat2, iat1) / lnlm.at(mu, l1);

        int sign = (mu % 2 == 0) ? 1 : -1;
        Complexf gam = static_cast<float>(sign) *
                       xclm_at(mn, l2, iat2, iat1) * lnlm.at(mu, l2);

        sum += gamtl * gam;
    }

    return sum;
}

// ============================================================================
// fms — main full multiple scattering routine
// ============================================================================
void fms(int lfms, int nsp, int ispin, int inclus, int npot,
         const Complexf* ck, const int* lipotx, const Complexf* xphase,
         int ik, int iverb, int minv, float rdirec,
         float toler1, float toler2, const bool* lcalc,
         std::vector<Eigen::MatrixXcf>& gg,
         FMSData& data) {

    const auto& cluster = data.cluster;
    const auto& rot     = data.rotation;
    const auto& lnlm    = data.lnlm;
    const auto& dw      = data.dw;
    auto& basis         = data.basis;
    const auto& cg      = data.cg;

    constexpr int nsp_lx2 = nspx * (lx + 1) * (lx + 1);

    // Helper to index lipotx (0-based potential index)
    std::vector<int> lipot(nphx + 1);
    for (int i = 0; i <= nphx; ++i) {
        lipot[i] = lipotx[i];
        if (lipot[i] <= 0) lipot[i] = lx;
        if (lipot[i] > lx)  lipot[i] = lx;
    }

    // Helper to access xphase(isp, l, iph) — l in [-lx, lx]
    // Fortran layout: xphase(nspx, -lx:lx, 0:nphasx)
    auto phase = [&](int isp, int l, int iph) -> Complexf {
        return xphase[isp + nspx * ((l + lx) + (2 * lx + 1) * iph)];
    };

    // Initialize gg output
    gg.resize(nphasx + 1);
    for (int ip = 0; ip <= nphasx; ++ip) {
        gg[ip] = Eigen::MatrixXcf::Zero(nsp_lx2, nsp_lx2);
    }

    // Determine potential range
    int ipi, ipf;
    if (lfms == 0) {
        // Extended system: only compute for one potential
        ipi = cluster.iphx[0];
        ipf = cluster.iphx[0];
    } else {
        ipi = 0;
        ipf = npot;
    }

    // Get basis kets
    std::vector<int> i0(nphx + 1, -1);
    getkts(nsp, inclus, lipot.data(), i0.data(), data);
    int istate = basis.istate;

    // Sanity check i0
    for (int ip = ipi; ip <= ipf; ++ip) {
        if (i0[ip] < 0) {
            throw std::runtime_error(
                "fms: Cannot find representative atom for potential " +
                std::to_string(ip) + ". Increase FMS radius.");
        }
    }

    // Runtime message
    if (iverb > 0 && minv == 0) {
        std::cout << "  FMS matrix (LUD) at point " << ik
                  << ", number of state kets = " << istate << std::endl;
    }

    // --- Compute c_lm(z) for all atom pairs ---
    // xclm[m][l][j][i] stored flat: (lx+1) x (lx+1) x nclusx x nclusx x nspx
    int xclm_size = (lx + 1) * (lx + 1) * nclusx * nclusx * nspx;
    std::vector<Complexf> xclm(xclm_size, Complexf(0.0f, 0.0f));

    auto xclm_idx = [](int m, int l, int j, int i, int isp) -> int {
        return m + (lx + 1) * (l + (lx + 1) * (j + nclusx * (i + nclusx * isp)));
    };

    // xrho[i][j][isp]
    std::vector<Complexf> xrho(nclusx * nclusx * nspx, Complexf(0.0f, 0.0f));
    auto rho_idx = [](int i, int j, int isp) -> int {
        return i + nclusx * (j + nclusx * isp);
    };

    // Temporary clm buffer
    constexpr int ltotb = lx + 1;
    constexpr int mntot = 2 * ltotb;
    constexpr int clm_rows = ltotb + 1;
    constexpr int clm_cols = mntot + 1;
    std::vector<Complexf> clm_buf(clm_rows * clm_cols);

    int lplus1 = lx + 1;
    int mplus1 = lx + 1;

    for (int i = 0; i < inclus; ++i) {
        for (int j = 0; j <= i; ++j) {
            // Distance between atoms i and j
            float r = 0.0f;
            for (int ix = 0; ix < 3; ++ix) {
                float d = cluster.rat(i, ix) - cluster.rat(j, ix);
                r += d * d;
            }
            r = std::sqrt(r);

            for (int isp = 0; isp < nsp; ++isp) {
                xrho[rho_idx(i, j, isp)] = ck[isp] * r;
                xrho[rho_idx(j, i, isp)] = xrho[rho_idx(i, j, isp)];

                if (i != j) {
                    xclmz(lplus1, mplus1, xrho[rho_idx(i, j, isp)], clm_buf.data());
                }
                for (int ll = 0; ll <= lx; ++ll) {
                    for (int mm = 0; mm <= lx; ++mm) {
                        if (i == j) {
                            xclm[xclm_idx(mm, ll, j, i, isp)] = Complexf(0.0f, 0.0f);
                        } else {
                            // Fortran: clm(ll+1, mm+1) -> C++: clm_buf[ll * clm_cols + mm]
                            Complexf val = clm_buf[ll * clm_cols + mm];
                            xclm[xclm_idx(mm, ll, j, i, isp)] = val;
                            xclm[xclm_idx(mm, ll, i, j, isp)] = val;
                        }
                    }
                }
            }
        }
    }

    // --- Fill G0 and T matrices ---
    Eigen::MatrixXcf g0  = Eigen::MatrixXcf::Zero(istate, istate);
    Eigen::MatrixXcf g0t = Eigen::MatrixXcf::Zero(istate, istate);
    Eigen::VectorXcf tmatrx_diag(istate);
    Eigen::VectorXcf tmatrx_offdiag(istate);
    tmatrx_diag.setZero();
    tmatrx_offdiag.setZero();

    Complexf coni_f(0.0f, 1.0f);
    float rdir2 = rdirec * rdirec;

    for (int ist1 = 0; ist1 < istate; ++ist1) {
        int iat1 = basis.iat(ist1);
        int l1   = basis.l(ist1);
        int m1   = basis.m(ist1);
        int isp1 = basis.isp(ist1);

        for (int ist2 = 0; ist2 < istate; ++ist2) {
            int iat2 = basis.iat(ist2);
            int l2   = basis.l(ist2);
            int m2   = basis.m(ist2);
            int isp2 = basis.isp(ist2);

            float rr = 0.0f;
            for (int ix = 0; ix < 3; ++ix) {
                float d = cluster.rat(iat1, ix) - cluster.rat(iat2, ix);
                rr += d * d;
            }

            if (iat1 == iat2) {
                // Same atom: G=0, compute T-matrix
                g0(ist1, ist2) = Complexf(0.0f, 0.0f);
                int iph = cluster.iphx[iat1];

                if (nsp == 1 && ispin == 0) {
                    if (ist1 == ist2) {
                        tmatrx_diag(ist1) =
                            (std::exp(2.0f * coni_f * phase(isp1, l1, iph)) - 1.0f)
                            / (2.0f * coni_f);
                    }
                } else {
                    if (ist1 == ist2) {
                        // Set spin index
                        int is = isp1;
                        if (nsp == 1) {
                            is = 0;
                            if (ispin > 0) is = 1;
                        }
                        // Diagonal T-matrix element
                        tmatrx_diag(ist1) =
                            (std::exp(2.0f * coni_f * phase(isp1, l1, iph)) - 1.0f)
                            / (2.0f * coni_f) * cg.jm(l1, m1, is) * cg.jm(l1, m1, is)
                          + (std::exp(2.0f * coni_f * phase(isp1, -l1, iph)) - 1.0f)
                            / (2.0f * coni_f) * cg.jp(l1, m1, is) * cg.jp(l1, m1, is);
                    } else if (nsp == 2 && l1 == l2 && m1 + isp1 == m2 + isp2) {
                        // Off-diagonal (spin-flip) T-matrix
                        // Fortran: tmatrx(nsp, ist1)
                        tmatrx_offdiag(ist1) =
                            (std::exp(2.0f * coni_f * phase(isp1, l1, iph)) - 1.0f
                           + std::exp(2.0f * coni_f * phase(isp2, l1, iph)) - 1.0f)
                            / (4.0f * coni_f)
                            * cg.jm(l1, m1, isp1) * cg.jm(l1, m2, isp2)
                          + (std::exp(2.0f * coni_f * phase(isp1, -l1, iph)) - 1.0f
                           + std::exp(2.0f * coni_f * phase(isp2, -l1, iph)) - 1.0f)
                            / (4.0f * coni_f)
                            * cg.jp(l1, m1, isp1) * cg.jp(l1, m2, isp2);
                    }
                }
            } else if (isp1 == isp2 && rr <= rdir2) {
                // Different atoms, same spin: compute G propagator
                g0(ist1, ist2) = Complexf(0.0f, 0.0f);

                // Pointer to xclm for this spin channel
                const Complexf* xclm_sp = &xclm[xclm_idx(0, 0, 0, 0, isp1)];

                for (int mu = -l1; mu <= l1; ++mu) {
                    int muabs = std::abs(mu);
                    Complexf gllmz = xgllm(muabs, ist1, ist2, xclm_sp, data);

                    // drix(mu, m1, l1, 1, iat2, iat1) * gllmz *
                    // drix(m2, mu, l2, 0, iat2, iat1)
                    g0(ist1, ist2) +=
                        rot.at(mu, m1, l1, 1, iat2, iat1) * gllmz *
                        rot.at(m2, mu, l2, 0, iat2, iat1);
                }

                // Spherical wave prefactor: exp(i*rho)/rho
                Complexf prefac = std::exp(coni_f * xrho[rho_idx(iat1, iat2, isp1)])
                                / xrho[rho_idx(iat1, iat2, isp1)];

                // Debye-Waller factor: exp(-sigma^2 * k^2 / bohr^2)
                // Fortran uses ck**2 (complex square), not |ck|^2 (modulus squared)
                float sig2_val = dw.at(iat1, iat2);
                Complexf ck_sq = ck[isp1] * ck[isp1];  // complex square, not std::norm
                prefac *= std::exp(-sig2_val * ck_sq
                                   / (static_cast<float>(feff::bohr) * static_cast<float>(feff::bohr)));

                g0(ist1, ist2) *= prefac;
            } else {
                // Different atoms, different spins: both G and T are zero
                g0(ist1, ist2) = Complexf(0.0f, 0.0f);
            }
        }
    }

    // --- Dispatch to solver ---
    int msord = 0;

    if (minv == 0) {
        gglu(nsp, i0.data(), ipi, ipf, lipot.data(), g0, tmatrx_diag,
             tmatrx_offdiag, g0t, gg, data);
    } else if (minv == 1) {
        ggbi(nsp, i0.data(), ipi, ipf, lipot.data(), g0, tmatrx_diag,
             tmatrx_offdiag, g0t, gg, toler1, toler2, lcalc, msord, data);
    } else if (minv == 2) {
        ggrm(nsp, i0.data(), ipi, ipf, lipot.data(), g0, tmatrx_diag,
             tmatrx_offdiag, g0t, gg, toler1, toler2, lcalc, msord, data);
    } else if (minv == 3) {
        gggm(nsp, i0.data(), ipi, ipf, lipot.data(), g0, tmatrx_diag,
             tmatrx_offdiag, g0t, gg, toler1, toler2, lcalc, msord, data);
    } else {
        ggtf(nsp, i0.data(), ipi, ipf, lipot.data(), g0, tmatrx_diag,
             tmatrx_offdiag, g0t, gg, toler1, toler2, lcalc, msord, data);
    }

    if (minv != 0) {
        std::string dec;
        switch (minv) {
            case 1: dec = "VdV"; break;
            case 2: dec = "LLU"; break;
            case 3: dec = "GMS"; break;
            default: dec = "TF"; break;
        }
        std::cout << "  Iterative FMS (" << dec << ") at point " << ik
                  << "; matrix size = " << istate
                  << "; MS order = " << msord << std::endl;
    }
}

} // namespace feff::fms
