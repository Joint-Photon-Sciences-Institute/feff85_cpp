// Radial integral routines for the ATOM module.
// Converted from: yzkteg.f, yzkrdf.f, dsordf.f, fdrirk.f, fdmocc.f, akeato.f

#include "radial_integrals.hpp"
#include "dsordf.hpp"
#include "coulomb_integrals.hpp"
#include "utility.hpp"
#include <cmath>
#include <algorithm>

namespace feff::atom {

// =========================================================================
// yzkteg — core yk/zk integration
// All arrays 0-based. Fortran indices shifted by -1.
// =========================================================================
void yzkteg(double f[], double af[], double g[], double ag[], const double dr[],
            double& ap, double h, int k, int nd, int np, int idim) {
    // Fortran: np = min(np, idim-2)
    np = std::min(np, idim - 2);

    double b = ap;
    ap = 0.0;
    g[0] = 0.0;
    g[1] = 0.0;

    // Development coefficients of yk at origin
    // Fortran: do 15 i=1,nd
    for (int i = 0; i < nd; ++i) {
        b += 1.0;
        ag[i] = af[i] / (b + k);
        if (af[i] != 0.0) {
            double c = std::pow(dr[0], b);
            g[0] += ag[i] * c;
            g[1] += ag[i] * std::pow(dr[1], b);
            af[i] = (k + k + 1) * ag[i] / (b - k - 1);
            ap += af[i] * c;
        }
    }

    // Multiply f by r
    for (int i = 0; i < np; ++i) {
        f[i] = f[i] * dr[i];
    }
    // Fortran: np1 = np+1 (1-based); C++: np is already 0-based count
    // f(np1) and f(np1+1) = 0 => f[np] and f[np+1] = 0
    f[np] = 0.0;
    f[np + 1] = 0.0;

    // Calculation of zk
    double eh = std::exp(h);
    double e = std::pow(eh, -k);
    b = h / 24.0;
    double c = 13.0 * b;
    double ee = e * e * b;
    b = b / e;

    // Fortran: do i=3,np1 => C++ i from 2 to np
    for (int i = 2; i <= np; ++i) {
        g[i] = g[i - 1] * e + (c * (f[i] + f[i - 1] * e) - (f[i - 2] * ee + f[i + 1] * b));
    }

    // Calculation of yk
    // Fortran: f(np) = g(np)
    f[np - 1] = g[np - 1];
    // Fortran: do i=np1,idim => f(i) = f(i-1)*e
    for (int i = np; i < idim; ++i) {
        f[i] = f[i - 1] * e;
    }

    // Fortran: i = k+k+1 (reuses variable i as integer constant)
    int kk1 = k + k + 1;
    b = kk1 * b * eh;
    ee = kk1 * ee / (eh * eh);
    e = e / eh;
    c = kk1 * c;

    // Fortran: do i=np-1,2,-1
    for (int i = np - 2; i >= 1; --i) {
        f[i] = f[i + 1] * e + (c * (g[i] + g[i + 1] * e) - (g[i + 2] * ee + g[i - 1] * b));
    }

    ee = e * e;
    c = 8.0 * c / 13.0;
    f[0] = f[2] * ee + c * (g[2] * ee + 4.0 * e * g[1] + g[0]);
    ap = (ap + f[0]) / std::pow(dr[0], k + 1);
}

// =========================================================================
// yzkrdf — calculate yk function for orbital products
// i, j are 0-based orbital indices (Fortran 1-based).
// When i < 0 (Fortran i <= 0), the function f is pre-constructed in work.dg.
// =========================================================================
void yzkrdf(int i, int j, int k,
            OrbitalArraysReal& orb, DiracWorkspaceReal& work,
            OrbitalConfig& config, InelasticFlag& inel, MeshParamsReal& mesh) {
    // NOTE: i and j follow Fortran convention (1-based orbital indices).
    // When i <= 0, j is a grid point count (not an orbital index).
    // When i > 0, both i and j are 1-based orbital indices.
    double chg[max_dev];

    if (i <= 0) {
        // Pre-constructed f in work.dg; dev coefficients in work.ag; power = k+2
        work.ap[0] = k + 2;
        // id = j directly (j is the integration limit, not an orbital index)
        int id = j;

        // Call yzkteg with work arrays
        yzkteg(work.dg, work.ag, work.dp, chg, mesh.dr,
               work.ap[0], mesh.hx, k, mesh.ndor, id, mesh.idim);
    } else {
        // Convert 1-based orbital indices to 0-based for C++ array access
        int i0 = i - 1;
        int j0 = j - 1;

        // Construct f(r) from orbital products
        double bgi[max_dev], bgj[max_dev], bpi[max_dev], bpj[max_dev];
        for (int l = 0; l < orb.ibgp; ++l) {
            bgi[l] = orb.bg[l][i0];
            bgj[l] = orb.bg[l][j0];
            bpi[l] = orb.bp[l][i0];
            bpj[l] = orb.bp[l][j0];
        }

        int id = std::min(config.nmax[i0], config.nmax[j0]);
        work.ap[0] = orb.fl[i0] + orb.fl[j0];

        if (inel.nem == 0) {
            // Full density: f = cg_i*cg_j + cp_i*cp_j
            for (int l = 0; l < id; ++l) {
                work.dg[l] = orb.cg[l][i0] * orb.cg[l][j0] + orb.cp[l][i0] * orb.cp[l][j0];
            }
            for (int l = 0; l < mesh.ndor; ++l) {
                work.ag[l] = aprdev(bgi, bgj, l + 1) + aprdev(bpi, bpj, l + 1);
            }
        } else {
            // Exchange: f = cg_i*cp_j
            for (int l = 0; l < id; ++l) {
                work.dg[l] = orb.cg[l][i0] * orb.cp[l][j0];
            }
            for (int l = 0; l < mesh.ndor; ++l) {
                work.ag[l] = aprdev(bgi, bpj, l + 1);
            }
        }

        yzkteg(work.dg, work.ag, work.dp, chg, mesh.dr,
               work.ap[0], mesh.hx, k, mesh.ndor, id, mesh.idim);
    }
}

// =========================================================================
// dsordf — Simpson overlap integral
// i, j are 0-based orbital indices.
// jnd selects the construction mode for the integrand.
// a = power of origin behavior (input for jnd >= 4 or negative jnd).
// =========================================================================
double dsordf(int i, int j, int n, int jnd, double a,
              OrbitalArraysReal& orb, DiracWorkspaceReal& work,
              OrbitalConfig& config, MeshParamsReal& mesh) {
    double hg[atom_grid] = {};
    double chg[max_dev] = {};

    double b = 0.0;
    int max0;
    bool do_construct_101 = false;

    if (jnd > 3) {
        // For jnd > 3, j is a point count (not orbital index)
        max0 = j;
        b = a;
        do_construct_101 = true;
    } else {
        max0 = std::min(config.nmax[i], config.nmax[j]);
        double bgi[max_dev], bgj[max_dev], bpi[max_dev], bpj[max_dev];
        for (int l = 0; l < orb.ibgp; ++l) {
            bgi[l] = orb.bg[l][i];
            bgj[l] = orb.bg[l][j];
            bpi[l] = orb.bp[l][i];
            bpj[l] = orb.bp[l][j];
        }

        if (std::abs(jnd) <= 2) {
            if (std::abs(jnd) < 2) {
                // jnd = +/-1: hg = cg_i*cg_j + cp_i*cp_j
                for (int l = 0; l < max0; ++l) {
                    hg[l] = orb.cg[l][i] * orb.cg[l][j] + orb.cp[l][i] * orb.cp[l][j];
                }
                for (int l = 0; l < mesh.ndor; ++l) {
                    chg[l] = aprdev(bgi, bgj, l + 1) + aprdev(bpi, bpj, l + 1);
                }
            } else {
                // jnd = +/-2: hg = cg_i*cp_j
                for (int l = 0; l < max0; ++l) {
                    hg[l] = orb.cg[l][i] * orb.cp[l][j];
                }
                for (int l = 0; l < mesh.ndor; ++l) {
                    chg[l] = aprdev(bgi, bpj, l + 1);
                }
            }
            b = orb.fl[i] + orb.fl[j];

            if (jnd <= 0) {
                // Multiply by dg (for negative jnd)
                for (int l = 0; l < max0; ++l) {
                    hg[l] = hg[l] * work.dg[l];
                }
                // Save chg into ap temporarily for the convolution
                for (int l = 0; l < mesh.ndor; ++l) {
                    work.ap[l] = chg[l];
                }
                b = b + a;
                for (int l = 0; l < mesh.ndor; ++l) {
                    chg[l] = aprdev(work.ap, work.ag, l + 1);
                }
            }
            // Fall through to integration (goto 301)
        } else {
            // jnd == 3: construct at label 101
            do_construct_101 = true;
        }
    }

    if (do_construct_101) {
        if (jnd == 4) {
            // hg = dg^2 + dp^2
            for (int l = 0; l < max0; ++l) {
                hg[l] = work.dg[l] * work.dg[l] + work.dp[l] * work.dp[l];
            }
            b = b + b;
            for (int l = 0; l < mesh.ndor; ++l) {
                chg[l] = aprdev(work.ag, work.ag, l + 1) + aprdev(work.ap, work.ap, l + 1);
            }
        } else if (jnd == 3) {
            // jnd == 3: hg = dg*cg_i + dp*cp_j
            max0 = std::min(config.nmax[i], config.nmax[j]);
            for (int l = 0; l < max0; ++l) {
                hg[l] = work.dg[l] * orb.cg[l][i] + work.dp[l] * orb.cp[l][j];
            }
            b = a + orb.fl[i];
            double bgi_loc[max_dev], bpj_loc[max_dev];
            for (int l = 0; l < orb.ibgp; ++l) {
                bgi_loc[l] = orb.bg[l][i];
                bpj_loc[l] = orb.bp[l][j];
            }
            for (int l = 0; l < mesh.ndor; ++l) {
                chg[l] = aprdev(bgi_loc, work.ag, l + 1) + aprdev(bpj_loc, work.ap, l + 1);
            }
        }
        // jnd >= 5: hg and chg are assumed pre-constructed (left as zero-initialized)
    }

    // Integration of hg * r^(n+1) by Simpson's method (label 301)
    double result = 0.0;
    int io = n + 1;
    for (int l = 0; l < max0; ++l) {
        hg[l] = hg[l] * std::pow(mesh.dr[l], io);
    }
    // Fortran: do l=2,max0,2 => dsordf += hg(l) + hg(l) + hg(l+1)
    // C++ 0-based: Fortran l=2 => C++ l=1; step 2; up to max0-1
    for (int l = 1; l < max0; l += 2) {
        result += hg[l] + hg[l] + hg[l + 1];
    }
    result = mesh.hx * (result + result + hg[0] - hg[max0 - 1]) / 3.0;

    // Integral from 0 to dr[0] using development coefficients
    b = b + n;
    for (int l = 0; l < mesh.ndor; ++l) {
        b += 1.0;
        result += chg[l] * std::pow(mesh.dr[0], b) / b;
    }

    return result;
}

// =========================================================================
// fdrirk — radial integrals Rk
// All orbital indices are 0-based.
// =========================================================================
double fdrirk(int i, int j, int l, int m, int k,
              OrbitalArraysReal& orb, DiracWorkspaceReal& work,
              OrbitalConfig& config, InelasticFlag& inel, MeshParamsReal& mesh) {
    double result = 0.0;

    if (i >= 0 && j >= 0) {
        // Fortran: i > 0 and j > 0  (1-based); C++: i >= 0 and j >= 0 (0-based)
        // yzkrdf expects 1-based orbital indices, so convert
        yzkrdf(i + 1, j + 1, k, orb, work, config, inel, mesh);

        int nn = std::abs(config.kap[i]) + std::abs(config.kap[j]);
        nn = std::max(nn - k, 1);
        double a_val = static_cast<double>(k + 1);

        // Shift ag coefficients: hg[nn-1+n-1] = -ag[n-1] for n=1..ndor
        // Then ag = hg, ag[0] += ap[0]
        double hg_tmp[max_dev] = {};
        for (int n = 0; n < mesh.ndor; ++n) {
            // Fortran: nn starts at abs(kap(i))+abs(kap(j))-k (or 1), 1-based
            // hg(nn) = -ag(n) where nn increments each iteration
            int idx = nn - 1 + n;  // 0-based target
            if (idx < mesh.ndor) {
                hg_tmp[idx] = -work.ag[n];
            }
        }
        for (int n = 0; n < mesh.ndor; ++n) {
            work.ag[n] = hg_tmp[n];
        }
        work.ag[0] += work.ap[0];

        // Now compute the overlap integral
        if (l < 0 || m < 0) return result;
        int n_jnd = -1;
        if (inel.nem != 0) n_jnd = -2;
        result = dsordf(l, m, -1, n_jnd, a_val, orb, work, config, mesh);
    } else {
        // i <= 0 case: just compute overlap
        if (l < 0 || m < 0) return result;
        double a_val = static_cast<double>(k + 1);
        int n_jnd = -1;
        if (inel.nem != 0) n_jnd = -2;
        result = dsordf(l, m, -1, n_jnd, a_val, orb, work, config, mesh);
    }

    return result;
}

// =========================================================================
// fdmocc — occupation number product
// i, j are 0-based orbital indices.
// =========================================================================
double fdmocc(int i, int j, OrbitalConfig& config) {
    if (j == i) {
        double result = config.xnel[i] * (config.xnel[j] - 1.0);
        double a = 2.0 * std::abs(config.kap[i]);
        return result * a / (a - 1.0);
    } else {
        return config.xnel[i] * config.xnel[j];
    }
}

// =========================================================================
// akeato — direct Coulomb angular coefficient
// Fortran: afgk(min(i,j), max(i,j), k/2)
// C++ struct: afgk[i][j][k] where first two dims are orbital indices
// =========================================================================
double akeato(int i, int j, int k, AngularCoefficients& ang) {
    if (i <= j) {
        return ang.afgk[i][j][k / 2];
    } else {
        return ang.afgk[j][i][k / 2];
    }
}

// =========================================================================
// bkeato — exchange Coulomb angular coefficient
// Fortran: afgk(max(i,j), min(i,j), k/2)
// Returns 0 when i == j.
// =========================================================================
double bkeato(int i, int j, int k, AngularCoefficients& ang) {
    if (i < j) {
        return ang.afgk[j][i][k / 2];
    } else if (i > j) {
        return ang.afgk[i][j][k / 2];
    }
    return 0.0;
}

// =========================================================================
// Wrapper overloads matching dsordf.hpp / coulomb_integrals.hpp signatures
// These delegate to the multi-argument versions above, unpacking AtomState.
// =========================================================================

double dsordf(int i, int j, int n, int jnd, double a, AtomState& state) {
    return dsordf(i, j, n, jnd, a, state.orb, state.work, state.config, state.mesh);
}

double fdrirk(int i, int j, int l, int m, int k, AtomState& state) {
    return fdrirk(i, j, l, m, k, state.orb, state.work, state.config, state.inelastic, state.mesh);
}

double akeato(int i, int j, int k, const AngularCoefficients& angular) {
    // const_cast is safe: the non-const version only reads from angular
    return akeato(i, j, k, const_cast<AngularCoefficients&>(angular));
}

double bkeato(int i, int j, int k, const AngularCoefficients& angular) {
    return bkeato(i, j, k, const_cast<AngularCoefficients&>(angular));
}

} // namespace feff::atom
