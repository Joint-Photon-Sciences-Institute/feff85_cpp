// Electron potential — converted from src/ATOM/potrdf.f

#include "potrdf.hpp"
#include "radial_integrals.hpp"
#include "utility.hpp"
#include <cmath>

namespace feff::atom {

void potrdf(int ia, AtomState& state) {
    auto& orb = state.orb;
    auto& work = state.work;
    auto& scf = state.scf;
    auto& config = state.config;
    auto& lag = state.lagrange;
    auto& nuc = state.nuclear;
    auto& mesh = state.mesh;
    auto& ang = state.angular;

    int ndor = mesh.ndor;
    int idim = mesh.idim;
    int norb = scf.norb;
    int norbsc = scf.norbsc;

    // Initialize development coefficients and working arrays
    for (int i = 0; i < ndor; ++i) {
        work.cep[i] = 0.0;
        work.ceg[i] = 0.0;
        work.av[i] = nuc.anoy[i];
    }

    double at[atom_grid] = {};
    double bt[atom_grid] = {};
    for (int i = 0; i < idim; ++i) {
        work.ep[i] = 0.0;
        work.eg[i] = 0.0;
        work.dv[i] = 0.0;
    }

    // Coulomb terms
    int jia = 2 * std::abs(config.kap[ia]) - 1;
    int k = 0;

    while (true) {
        // Zero workspace for density accumulation
        for (int i = 0; i < idim; ++i) work.dg[i] = 0.0;
        for (int i = 0; i < ndor; ++i) work.ag[i] = 0.0;

        int max0 = 0;
        for (int j = 0; j < norb; ++j) {
            double bgj[max_dev], bpj[max_dev];
            for (int i = 0; i < 10; ++i) {
                bgj[i] = orb.bg[i][j];
                bpj[i] = orb.bp[i][j];
            }

            int m = 2 * std::abs(config.kap[j]) - 1;
            if (k <= m) {
                double a = akeato(ia, j, k, ang) / config.xnel[ia];
                if (a != 0.0) {
                    m = config.nmax[j];
                    for (int i = 0; i < m; ++i) {
                        work.dg[i] += a * (orb.cg[i][j] * orb.cg[i][j] +
                                           orb.cp[i][j] * orb.cp[i][j]);
                    }
                    int n = 2 * std::abs(config.kap[j]) - k;
                    int l = ndor + 2 - n;
                    if (l > 0) {
                        double a2 = a * orb.fix[j] * orb.fix[j];
                        for (int i = 0; i < l; ++i) {
                            int m2 = n - 2 + i;  // 0-based
                            if (m2 >= 0 && m2 < ndor) {
                                work.ag[m2] += a2 * (aprdev(bgj, bgj, i + 1) +
                                                     aprdev(bpj, bpj, i + 1));
                            }
                        }
                    }
                }
            }
            if (config.nmax[j] > max0) max0 = config.nmax[j];
        }

        yzkrdf(0, max0, k, orb, work, config, state.inelastic, mesh);

        for (int i = 0; i < ndor; ++i) {
            int l = k + i + 3;  // 0-based: Fortran l=k+i+3
            if (l < ndor) {
                work.av[l] -= work.ag[i];
            }
        }
        for (int i = 0; i < idim; ++i) {
            work.dv[i] += work.dg[i];
        }

        k += 2;
        if (k <= ndor) {
            if (k - 1 < ndor) work.av[k - 1] += work.ap[0];
        }
        if (k >= jia) break;
    }

    // Exchange terms
    if (mesh.method != 0) {
        for (int j = 0; j < norb; ++j) {
            if (j == ia) continue;
            int max0 = config.nmax[j];
            int jj = 2 * std::abs(config.kap[j]) - 1;
            int kma = (jj + jia) / 2;
            k = std::abs(jj - kma);
            if (config.kap[j] * config.kap[ia] < 0) k += 1;

            while (k <= kma) {
                double a = bkeato(ia, j, k, ang) / config.xnel[ia];
                if (a != 0.0) {
                    yzkrdf(j + 1, ia + 1, k, orb, work, config, state.inelastic, mesh);
                    for (int i = 0; i < max0; ++i) {
                        work.eg[i] += a * work.dg[i] * orb.cg[i][j];
                        work.ep[i] += a * work.dg[i] * orb.cp[i][j];
                    }

                    int n = k + 1 + std::abs(config.kap[j]) - std::abs(config.kap[ia]);
                    if (n > 0 && n <= ndor) {
                        for (int i = n - 1; i < ndor; ++i) {
                            double coeff = a * work.ap[0] * orb.fix[j] / orb.fix[ia];
                            work.ceg[i] += orb.bg[i + 1 - n][j] * coeff;
                            work.cep[i] += orb.bp[i + 1 - n][j] * coeff;
                        }
                    }

                    int ii = 2 * std::abs(config.kap[j]) + 1;
                    if (ii <= ndor) {
                        double bgj[max_dev], bpj[max_dev];
                        for (int ix = 0; ix < 10; ++ix) {
                            bgj[ix] = orb.bg[ix][j];
                            bpj[ix] = orb.bp[ix][j];
                        }
                        for (int n2 = ii - 1; n2 < ndor; ++n2) {
                            double fix2 = orb.fix[j] * orb.fix[j];
                            work.ceg[n2] -= a * aprdev(work.ag, bgj, n2 + 2 - ii) * fix2;
                            work.cep[n2] -= a * aprdev(work.ag, bpj, n2 + 2 - ii) * fix2;
                        }
                    }
                }
                k += 2;
            }
        }
    }

    // Lagrange multiplier terms
    if (lag.ipl != 0) {
        for (int j = 0; j < norbsc; ++j) {
            if (config.kap[j] != config.kap[ia] || j == ia) continue;
            if (lag.nre[j] < 0 && lag.nre[ia] < 0) continue;
            int m_idx = std::max(j, ia);
            int i_idx = std::min(j, ia) + (m_idx * (m_idx - 1)) / 2;
            double a = lag.eps[i_idx] * config.xnel[j];
            int max0 = config.nmax[j];
            for (int i = 0; i < max0; ++i) {
                at[i] += a * orb.cg[i][j];
                bt[i] += a * orb.cp[i][j];
            }
            for (int i = 0; i < ndor; ++i) {
                work.ceg[i] += orb.bg[i][j] * a;
                work.cep[i] += orb.bp[i][j] * a;
            }
        }
    }

    // Add nuclear potential and divide by speed of light
    for (int i = 0; i < ndor; ++i) {
        work.av[i] /= work.cl;
        work.cep[i] /= work.cl;
        work.ceg[i] /= work.cl;
    }
    for (int i = 0; i < idim; ++i) {
        work.dv[i] = (work.dv[i] / mesh.dr[i] + nuc.dvn[i]) / work.cl;
        work.ep[i] = (work.ep[i] + bt[i] * mesh.dr[i]) / work.cl;
        work.eg[i] = (work.eg[i] + at[i] * mesh.dr[i]) / work.cl;
    }
}

} // namespace feff::atom
