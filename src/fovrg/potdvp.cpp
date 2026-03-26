// Coulomb potential development coefficients for FOVRG.
// Converted from: src/FOVRG/potdvp.f

#include "potdvp.hpp"
#include "radial_integrals_c.hpp"
#include <feff/dimensions.hpp>
#include <cmath>

namespace feff::fovrg {

static inline double fpow_int(double x, int n) {
    if (n == 0) return 1.0;
    if (n < 0) { x = 1.0/x; n = -n; }
    double r = 1.0;
    for (int i = 0; i < n; ++i) r *= x;
    return r;
}

void potdvp(FovrgState& state)
{
    auto& orb = state.orb;
    auto& config = state.config;
    auto& scf = state.scf;
    auto& work = state.work;
    auto& nuclear = state.nuclear;
    auto& mesh = state.mesh;
    int norb = scf.norb;
    int ndor = mesh.ndor;
    double cl = work.cl;

    // av = anoy (nuclear potential development coefficients)
    for (int i = 0; i < 10; i++) {
        work.av[i] = FeffComplex(nuclear.anoy[i], 0.0);
    }

    // Calculate density development coefficients
    // ag is used as workspace here
    for (int i = 0; i < ndor; i++) {
        work.ag[i] = FeffComplex(0.0, 0.0);
    }

    for (int j = 0; j < norb - 1; j++) {
        double bgj[10], bpj[10];
        for (int i = 0; i < 10; i++) {
            bgj[i] = orb.bg[i][j];
            bpj[i] = orb.bp[i][j];
        }
        int n = 2 * std::abs(config.kap[j]);
        int l = ndor + 2 - n;
        if (l > 0) {
            for (int i = 0; i < l; i++) {
                int m = n - 2 + i;  // 0-based index into ag
                work.ag[m] = work.ag[m] + config.xnel[j] *
                    FeffComplex(aprdep(bgj, bgj, i + 1) + aprdep(bpj, bpj, i + 1), 0.0) *
                    orb.fix[j] * orb.fix[j];
            }
        }
    }

    // Transform density coefficients into potential coefficients
    work.ap[0] = FeffComplex(0.0, 0.0);
    for (int i = 0; i < ndor; i++) {
        work.ag[i] = work.ag[i] / FeffComplex((i + 3) * (i + 2), 0.0);
        work.ap[0] = work.ap[0] + work.ag[i] * fpow_int(mesh.dr[0], i + 2);
    }

    for (int i = 0; i < ndor; i++) {
        int l = i + 3;  // 0-based target index in av
        if (l < ndor) {
            work.av[l] = work.av[l] - work.ag[i];
        }
    }

    // av(2) += ap(1) (add density integral at first point)
    work.av[1] = work.av[1] + work.ap[0];

    // Division by speed of light
    for (int i = 0; i < 10; i++) {
        work.av[i] = work.av[i] / cl;
    }
}

} // namespace feff::fovrg
