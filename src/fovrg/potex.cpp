// Exchange potential for photoelectron in FOVRG.
// Converted from: src/FOVRG/potex.f

#include "potex.hpp"
#include "radial_integrals_c.hpp"
#include <feff/dimensions.hpp>
#include <cmath>

namespace feff::fovrg {

void potex(const FeffComplex ps[], const FeffComplex qs[],
           const FeffComplex aps[10], const FeffComplex aqs[10],
           int jri, FovrgState& state)
{
    auto& orb = state.orb;
    auto& config = state.config;
    auto& scf = state.scf;
    auto& work = state.work;
    auto& mesh = state.mesh;
    auto& angular = state.angular;
    int norb = scf.norb;
    int ndor = mesh.ndor;
    int idim = mesh.idim;
    double cl = work.cl;

    // Zero exchange potentials
    for (int i = 0; i < 10; i++) {
        work.cep[i] = FeffComplex(0.0, 0.0);
        work.ceg[i] = FeffComplex(0.0, 0.0);
    }
    for (int i = 0; i < idim; i++) {
        work.ep[i] = FeffComplex(0.0, 0.0);
        work.eg[i] = FeffComplex(0.0, 0.0);
    }

    // jia = 2*|kap(norb)| - 1
    int jia = 2 * std::abs(config.kap[norb - 1]) - 1;

    // Exchange terms: loop over core orbitals
    for (int j = 0; j < norb - 1; j++) {
        int jj = 2 * std::abs(config.kap[j]) - 1;
        int kma = (jj + jia) / 2;
        int k = std::abs(jj - kma);
        if ((config.kap[j] * config.kap[norb - 1]) < 0) k = k + 1;
        int kmin = k;

        while (k <= kma) {
            double a = angular(config.kap[norb - 1], j, (k - kmin) / 2);
            if (a != 0.0) {
                // Calculate yk integral
                yzkrdc(j, k, orb.fl[norb - 1], ps, qs, aps, aqs, state);

                // Accumulate exchange potential
                // work.gg holds dg output from yzkrdc, work.gp holds dp
                for (int i = 0; i < idim; i++) {
                    work.eg[i] = work.eg[i] + a * work.gg[i] * orb.cg[i][j];
                    work.ep[i] = work.ep[i] + a * work.gg[i] * orb.cp[i][j];
                }

                // Development coefficients
                // n = k+1+|kap(j)|-|kap(norb)|
                int n = k + 1 + std::abs(config.kap[j]) - std::abs(config.kap[norb - 1]);
                // Different for irregular solution
                if (orb.fl[norb - 1] < 0.0)
                    n = k + 1 + std::abs(config.kap[j]) + std::abs(config.kap[norb - 1]);

                if (n <= ndor) {
                    for (int i = n - 1; i < ndor; i++) {
                        // Fortran: ceg(i) += bg(i+1-n,j)*a*ap(1)*fix(j)/fix(norb)
                        // 0-based: bg[i-n+1][j], ap[0]
                        work.ceg[i] = work.ceg[i] + orb.bg[i - n + 1][j] * a *
                                      work.ap[0] * orb.fix[j] / orb.fix[norb - 1];
                        work.cep[i] = work.cep[i] + orb.bp[i - n + 1][j] * a *
                                      work.ap[0] * orb.fix[j] / orb.fix[norb - 1];
                    }
                }

                // Second correction term
                int i_start = 2 * std::abs(config.kap[j]) + 1;
                if (i_start <= ndor) {
                    double bgj[10], bpj[10];
                    for (int ix = 0; ix < 10; ix++) {
                        bgj[ix] = orb.bg[ix][j];
                        bpj[ix] = orb.bp[ix][j];
                    }
                    for (int nn = i_start - 1; nn < ndor; nn++) {
                        // Fortran: nx = n + 1 - i where n is loop var, i is i_start
                        // 0-based: nx = nn - i_start + 2 => 1-based: nx = nn+1 - i_start + 1
                        int nx = nn - i_start + 2;
                        work.ceg[nn] = work.ceg[nn] - a * aprdec(work.ag, bgj, nx) *
                                       orb.fix[j] * orb.fix[j];
                        work.cep[nn] = work.cep[nn] - a * aprdec(work.ag, bpj, nx) *
                                       orb.fix[j] * orb.fix[j];
                    }
                }
            }
            k = k + 2;
        }
    }

    // Division by speed of light
    for (int i = 0; i < ndor; i++) {
        work.cep[i] = work.cep[i] / cl;
        work.ceg[i] = work.ceg[i] / cl;
    }
    for (int i = 0; i < jri; i++) {
        work.ep[i] = work.ep[i] / cl;
        work.eg[i] = work.eg[i] / cl;
    }
    for (int i = jri; i < nrptx; i++) {
        work.ep[i] = FeffComplex(0.0, 0.0);
        work.eg[i] = FeffComplex(0.0, 0.0);
    }
}

} // namespace feff::fovrg
