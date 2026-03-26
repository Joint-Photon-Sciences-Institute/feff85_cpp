// Geometric preparation for FMS calculation.
// Converted from: xprep.f and helper routines in xstaff.f
// Original authors: Bruce Ravel, Alexei Ankudinov
//
// Key changes from Fortran:
//   - COMMON blocks replaced by FMSData struct references
//   - 1-based indexing converted to 0-based
//   - Debye-Waller model calls (sigms, sigem, sigrm) are stubs
//     pending conversion of the DEBYE module

#include "xprep.hpp"

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <cstring>

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>

#include "fms_types.hpp"
#include "../math/wigner.hpp"

namespace feff::fms {

// ============================================================================
// xfctst — factorial table for Legendre normalization
// Adapted from xstaff.f xfctst
// ============================================================================
void xfctst(float& afac, float& flzero, float* flg) {
    afac = 0.03125f;
    flzero = 1.0f;
    flg[0] = 1.0f;
    flg[1] = afac;
    for (int i = 2; i <= 50; ++i) {
        flg[i] = flg[i - 1] * i * afac;
    }
}

// ============================================================================
// xanlm — Legendre polynomial normalization factors
// xnlm(m,l) = sqrt((2l+1)(l-m)!/(l+m)!) * afac^m
// ============================================================================
void xanlm(int lmaxp1, int mmaxp1, LegendreNorm& lnlm) {
    float afac, flzero;
    float flg[51];
    xfctst(afac, flzero, flg);

    for (int il = 1; il <= lmaxp1; ++il) {
        int mmxp1 = std::min(mmaxp1, il);
        for (int im = 1; im <= mmxp1; ++im) {
            int l = il - 1;
            int m = im - 1;
            float cnlm = (2 * l + 1) * flg[l - m] / flg[l + m];
            cnlm = std::sqrt(cnlm) * std::pow(afac, m);
            lnlm.at(m, l) = cnlm;
        }
    }
}

// ============================================================================
// atheap — heap-sort atoms by distance from central atom (origin)
// Adapted from xstaff.f atheap.  rat is stored as rat[iat*3 + dim].
// ============================================================================
void atheap(int nat, float* rat, int* iphat, double* ra) {
    if (nat < 2) return;

    // Compute squared distances with small perturbation to stabilize sort order
    int needs_sort = 0;
    for (int i = 0; i < nat; ++i) {
        double r2 = static_cast<double>(rat[i * 3 + 0]) * rat[i * 3 + 0]
                   + static_cast<double>(rat[i * 3 + 1]) * rat[i * 3 + 1]
                   + static_cast<double>(rat[i * 3 + 2]) * rat[i * 3 + 2];
        ra[i] = r2 + (i + 1) * 1.0e-8;
        if (needs_sort == 0 && i > 0 && ra[i] < ra[i - 1]) {
            needs_sort = 1;
        }
    }
    if (!needs_sort) return;

    // Heapsort (adapted from Numerical Recipes)
    int l = nat / 2;
    int ir = nat - 1;

    float toss[3];
    int itoss;
    double dum;

    while (true) {
        if (l > 0) {
            --l;
            toss[0] = rat[l * 3 + 0];
            toss[1] = rat[l * 3 + 1];
            toss[2] = rat[l * 3 + 2];
            itoss = iphat[l];
            dum = ra[l];
        } else {
            toss[0] = rat[ir * 3 + 0];
            toss[1] = rat[ir * 3 + 1];
            toss[2] = rat[ir * 3 + 2];
            itoss = iphat[ir];
            dum = ra[ir];

            rat[ir * 3 + 0] = rat[0];
            rat[ir * 3 + 1] = rat[1];
            rat[ir * 3 + 2] = rat[2];
            iphat[ir] = iphat[0];
            ra[ir] = ra[0];

            --ir;
            if (ir == 0) {
                rat[0] = toss[0];
                rat[1] = toss[1];
                rat[2] = toss[2];
                iphat[0] = itoss;
                ra[0] = dum;
                return;
            }
        }

        int i = l;
        int j = l + l + 1;  // 0-based: left child

        while (j <= ir) {
            if (j < ir && ra[j] < ra[j + 1]) {
                ++j;
            }
            if (dum < ra[j]) {
                rat[i * 3 + 0] = rat[j * 3 + 0];
                rat[i * 3 + 1] = rat[j * 3 + 1];
                rat[i * 3 + 2] = rat[j * 3 + 2];
                iphat[i] = iphat[j];
                ra[i] = ra[j];
                i = j;
                j = 2 * j + 1;
            } else {
                j = ir + 1;
            }
        }

        rat[i * 3 + 0] = toss[0];
        rat[i * 3 + 1] = toss[1];
        rat[i * 3 + 2] = toss[2];
        iphat[i] = itoss;
        ra[i] = dum;
    }
}

// ============================================================================
// getang — polar angles of vector R_i - R_j
// ============================================================================
void getang(const float* xrat, int nclusx_dim, int i, int j,
            float& theta, float& phi_out) {
    constexpr float tiny = 1.0e-7f;
    constexpr float pi_f = static_cast<float>(feff::pi);

    float x = xrat[i * 3 + 0] - xrat[j * 3 + 0];
    float y = xrat[i * 3 + 1] - xrat[j * 3 + 1];
    float z = xrat[i * 3 + 2] - xrat[j * 3 + 2];
    float r = std::sqrt(x * x + y * y + z * z);

    phi_out = 0.0f;
    theta = 0.0f;

    if (i != j) {
        if (std::abs(x) < tiny) {
            if (std::abs(y) < tiny)
                phi_out = 0.0f;
            else if (y > tiny)
                phi_out = pi_f / 2.0f;
            else
                phi_out = -pi_f / 2.0f;
        } else {
            phi_out = std::atan2(y, x);
        }

        if (r > tiny) {
            if (z <= -r)
                theta = pi_f;
            else if (z < r)
                theta = std::acos(z / r);
            // else theta = 0 (already set)
        }
    }
}

// ============================================================================
// rotint — initialize rotation matrix save cache
// ============================================================================
void rotint(RotationData& rot) {
    rot.jsav = 0;
    for (int js = 0; js < RotationData::jsavx; ++js) {
        rot.betsav[js] = static_cast<float>(RotationData::jbmagk);
        rot.ldsav[js] = 0;
        rot.mdsav[js] = 0;
    }
    std::fill(rot.drisav.begin(), rot.drisav.end(), 0.0f);
}

// ============================================================================
// rotxan — compute rotation matrix for atom pair (i,j)
// Adapted from xstaff.f rotxan.  Edmonds' recursion.
// ============================================================================
void rotxan(int lxp1, int mxp1, float betax, int i, int j, int k,
            ClusterData& cluster, RotationData& rot) {
    constexpr float pi_f = static_cast<float>(feff::pi);
    constexpr int lxx = 24;

    // Check cache for previously computed rotation matrix
    for (int isav = 0; isav < rot.jsav; ++isav) {
        if (rot.betsav[isav] == RotationData::jbmagk) break;
        if (lxp1 == rot.ldsav[isav] && mxp1 == rot.mdsav[isav] &&
            std::abs(betax - rot.betsav[isav]) <= RotationData::roteps) {
            // Use cached value
            for (int il = 0; il <= lx; ++il) {
                for (int m1 = -il; m1 <= il; ++m1) {
                    for (int m2 = -il; m2 <= il; ++m2) {
                        rot.at(m2, m1, il, k, j, i) =
                            Complexf(rot.drisav[rot.drisav_index(m2, m1, il, isav)], 0.0f);
                    }
                }
            }
            goto apply_phi;
        }
    }

    {
        // Compute rotation matrix from scratch using Edmonds recursion
        // dri0 indexed as dri0[il][m][n], il in [0,lxx], m,n in [0, 2*lxx]
        constexpr int dri_dim = 2 * lxx + 1;
        std::vector<float> dri0((lxx + 1) * dri_dim * dri_dim, 0.0f);

        auto dri = [&](int il, int m, int n) -> float& {
            return dri0[il + (lxx + 1) * (m + dri_dim * n)];
        };

        int nm = mxp1;
        int ndm = lxp1 + nm - 1;
        float xc = std::cos(betax / 2.0f);
        float xs = std::sin(betax / 2.0f);
        float s  = std::sin(betax);

        dri(0, 0, 0) = 1.0f;
        dri(1, 0, 0) = xc * xc;
        dri(1, 0, 1) = s / std::sqrt(2.0f);
        dri(1, 0, 2) = xs * xs;
        dri(1, 1, 0) = -dri(1, 0, 1);
        dri(1, 1, 1) = std::cos(betax);
        dri(1, 1, 2) = dri(1, 0, 1);
        dri(1, 2, 0) = dri(1, 0, 2);
        dri(1, 2, 1) = -dri(1, 1, 2);
        dri(1, 2, 2) = dri(1, 0, 0);

        for (int il = 2; il < lxp1; ++il) {
            int ln = std::min(2 * il + 1, ndm);
            int lm = std::min(2 * il - 1, ndm);
            for (int n = 0; n < ln; ++n) {
                for (int m = 0; m < lm; ++m) {
                    float t1 = static_cast<float>((2 * il - n) * (2 * il - 1 - n));
                    float t  = static_cast<float>((2 * il - m) * (2 * il - 1 - m));
                    float f1 = std::sqrt(t1 / t);
                    float f2 = std::sqrt(static_cast<float>((2 * il - n) * n) / t);
                    float t3 = static_cast<float>((n - 1) * n);
                    float f3 = std::sqrt(t3 / t);
                    float dlnm = f1 * xc * xc * dri(il - 1, n, m);
                    if (n > 0) dlnm -= f2 * s * dri(il - 1, n - 1, m);
                    if (n > 1) dlnm += f3 * xs * xs * dri(il - 1, n - 2, m);
                    dri(il, n, m) = dlnm;
                    if (n >= lm) {
                        int sign = ((n - m) % 2 == 0) ? 1 : -1;
                        dri(il, m, n) = sign * dri(il, n, m);
                    }
                }
                if (n >= lm) {
                    int ll = 2 * il;
                    dri(il, ll - 1, ll - 1) = dri(il, 1, 1);
                    dri(il, ll, ll - 1) = -dri(il, 0, 1);
                    dri(il, ll - 1, ll) = -dri(il, 1, 0);
                    dri(il, ll, ll) = dri(il, 0, 0);
                }
            }
        }

        // Initialize drix for this pair
        for (int il = 0; il <= lx; ++il) {
            for (int m1 = -lx; m1 <= lx; ++m1) {
                for (int m2 = -lx; m2 <= lx; ++m2) {
                    rot.at(m2, m1, il, k, j, i) = Complexf(0.0f, 0.0f);
                    rot.at(m2, m1, il, k, i, i) = Complexf(0.0f, 0.0f);
                }
            }
        }

        // Copy dri0 into drix
        for (int il = 0; il < lxp1; ++il) {
            int mmx = std::min(il, mxp1 - 1);
            for (int m1 = -mmx; m1 <= mmx; ++m1) {
                for (int m2 = -mmx; m2 <= mmx; ++m2) {
                    rot.at(m2, m1, il, k, j, i) =
                        Complexf(dri(il, m1 + il, m2 + il), 0.0f);
                }
            }
        }

        // Save to cache if room
        if (rot.jsav < RotationData::jsavx) {
            int idx = rot.jsav++;
            rot.betsav[idx] = betax;
            rot.ldsav[idx] = lxp1;
            rot.mdsav[idx] = mxp1;
            for (int il = 0; il <= lx; ++il) {
                for (int m1 = -il; m1 <= il; ++m1) {
                    for (int m2 = -il; m2 <= il; ++m2) {
                        rot.drisav[rot.drisav_index(m2, m1, il, idx)] =
                            rot.at(m2, m1, il, k, j, i).real();
                    }
                }
            }
        }
    }

apply_phi:
    // Apply azimuthal phase factor exp(i*m*(phi-pi))
    {
        Complexf coni_f(0.0f, 1.0f);
        float phi_ij = cluster.phi(i, j);

        for (int il = 0; il <= lx; ++il) {
            for (int m1 = -il; m1 <= il; ++m1) {
                Complexf dum = coni_f * static_cast<float>(m1) * (phi_ij - pi_f);
                if (k == 1) dum = -dum;
                dum = std::exp(dum);
                for (int m2 = -il; m2 <= il; ++m2) {
                    if (k == 1) {
                        rot.at(m2, m1, il, k, j, i) *= dum;
                    } else {
                        rot.at(m1, m2, il, k, j, i) *= dum;
                    }
                }
            }
        }
    }
}

// ============================================================================
// xprep — main geometric preparation
// ============================================================================
void xprep(int iph0, int idwopt, int nat, int& inclus, int npot,
           const int* iphat, float rmax, const float* rat,
           const int* izx, float rnrmav, float temper, float thetad, float sig2,
           int minv, float rdirec,
           FMSData& data) {

    auto& cluster = data.cluster;
    auto& rot     = data.rotation;
    auto& lnlm    = data.lnlm;
    auto& dw      = data.dw;
    auto& cg      = data.cg;

    // Initialize geometrical arrays
    std::fill(cluster.xphi.begin(), cluster.xphi.end(), 0.0f);
    std::fill(cluster.xrat.begin(), cluster.xrat.end(), 0.0f);
    std::fill(cluster.iphx.begin(), cluster.iphx.end(), 0);
    inclus = 0;

    // Working copies of rat and iphat (0-based)
    std::vector<float> rat2(nat * 3);
    std::vector<int> iphat2(nat);

    // Find the central atom (ipot == iph0)
    int icen = -1;
    for (int i = 0; i < nat; ++i) {
        iphat2[i] = iphat[i];
        if (iphat[i] == iph0) {
            if (icen < 0) {
                icen = i;
            } else if (iph0 == 0) {
                throw std::runtime_error(
                    "xprep: More than one atom with ipot=0 in extended cluster");
            }
        }
    }
    if (icen < 0) {
        throw std::runtime_error("xprep: No central atom found with ipot=" +
                                 std::to_string(iph0));
    }

    // Center coordinates on central atom
    for (int i = 0; i < nat; ++i) {
        rat2[i * 3 + 0] = rat[i * 3 + 0] - rat[icen * 3 + 0];
        rat2[i * 3 + 1] = rat[i * 3 + 1] - rat[icen * 3 + 1];
        rat2[i * 3 + 2] = rat[i * 3 + 2] - rat[icen * 3 + 2];
    }

    // Sort atoms by distance from central atom
    std::vector<double> ra(nat);
    atheap(nat, rat2.data(), iphat2.data(), ra.data());

    // Define cluster: atoms within rmax of center
    float rmax2 = rmax * rmax;
    inclus = 0;
    for (int i = 0; i < nat; ++i) {
        float rr = rat2[i * 3 + 0] * rat2[i * 3 + 0]
                 + rat2[i * 3 + 1] * rat2[i * 3 + 1]
                 + rat2[i * 3 + 2] * rat2[i * 3 + 2];
        if (rr > rmax2) {
            inclus = i;
            break;
        }
    }
    if (inclus == 0) inclus = nat;

    // Cap at nclusx
    if (inclus > nclusx) {
        inclus = nclusx;
    }

    // Copy into cluster structure
    for (int iat = 0; iat < inclus; ++iat) {
        cluster.iphx[iat] = iphat2[iat];
        cluster.rat(iat, 0) = rat2[iat * 3 + 0];
        cluster.rat(iat, 1) = rat2[iat * 3 + 1];
        cluster.rat(iat, 2) = rat2[iat * 3 + 2];
    }

    // Calculate rotation matrix elements and phi angles
    rotint(rot);
    int lplus1 = lx + 1;
    int mplus1 = lx + 1;

    for (int i = 0; i < inclus; ++i) {
        for (int j = 0; j < inclus; ++j) {
            float xbeta, xphi_val;
            getang(cluster.xrat.data(), nclusx, i, j, xbeta, xphi_val);
            cluster.phi(i, j) = xphi_val;

            if (i != j) {
                for (int kk = 0; kk <= 1; ++kk) {
                    float beta = (kk == 1) ? -xbeta : xbeta;
                    rotxan(lplus1, mplus1, beta, i, j, kk, cluster, rot);
                }
            }
        }
    }

    // Calculate spherical harmonic normalization factors
    xanlm(lplus1, mplus1, lnlm);

    // Initialize Debye-Waller factors
    std::fill(dw.sigsqr.begin(), dw.sigsqr.end(), 0.0f);

    if (idwopt >= 0) {
        // Calculate Debye-Waller factors for each atom pair
        // For CD model (idwopt==0), we would call sigms here.
        // Currently we just add the global sig2 contribution.
        for (int iat1 = 0; iat1 < inclus - 1; ++iat1) {
            for (int iat2 = iat1 + 1; iat2 < inclus; ++iat2) {
                // TODO: Call Debye-Waller model (sigms/sigem/sigrm) when
                // the DEBYE module is converted to C++.
                // For now, apply only the global sigma^2.
                dw.at(iat1, iat2) = sig2;
                dw.at(iat2, iat1) = sig2;
            }
        }
    }

    // Calculate Clebsch-Gordan coefficients <LS|J>
    for (int l1 = 0; l1 <= lx; ++l1) {
        for (int mm = -l1; mm <= l1; ++mm) {
            for (int isp1 = 0; isp1 < 2; ++isp1) {
                int j1 = 2 * l1;
                int j2 = 1;
                int j3p = j1 + 1;
                int j3m = j1 - 1;
                int m1 = 2 * mm;
                int m2 = 2 * (isp1 + 1) - 3;  // isp1=0 -> -1, isp1=1 -> +1

                // j = l + 1/2
                float val_p = std::sqrt(static_cast<float>(j3p + 1)) *
                    static_cast<float>(feff::math::cwig3j(j1, j2, j3p, m1, m2, 2));
                if (((j2 - j1 - m1 - m2) / 2) % 2 != 0)
                    val_p = -val_p;
                cg.jp(l1, mm, isp1) = val_p;

                // j = l - 1/2
                if (j3m >= 0) {
                    float val_m = std::sqrt(static_cast<float>(j3m + 1)) *
                        static_cast<float>(feff::math::cwig3j(j1, j2, j3m, m1, m2, 2));
                    if (((j2 - j1 - m1 - m2) / 2) % 2 != 0)
                        val_m = -val_m;
                    cg.jm(l1, mm, isp1) = val_m;
                }
            }
        }
    }
}

} // namespace feff::fms
