// Schmidt orthogonalization for FOVRG.
// Converted from: src/FOVRG/ortdac.f

#include "ortdac.hpp"
#include "radial_integrals_c.hpp"
#include <feff/dimensions.hpp>

namespace feff::fovrg {

void ortdac(int ikap, FeffComplex ps[], FeffComplex qs[],
            FeffComplex aps[10], FeffComplex aqs[10], FovrgState& state)
{
    auto& orb = state.orb;
    auto& config = state.config;
    auto& scf = state.scf;
    auto& mesh = state.mesh;
    int norb = scf.norb;
    int ndor = mesh.ndor;
    int idim = mesh.idim;

    for (int j = 0; j < norb - 1; j++) {
        if (config.kap[j] == ikap && config.xnel[j] > 0.0) {
            // Calculate overlap integral
            FeffComplex a = dsordc(j, FeffComplex(orb.fl[norb - 1], 0.0),
                                   ps, qs, aps, aqs, orb, mesh);

            // Subtract projection
            for (int i = 0; i < idim; i++) {
                ps[i] = ps[i] - a * orb.cg[i][j];
                qs[i] = qs[i] - a * orb.cp[i][j];
            }
            for (int i = 0; i < ndor; i++) {
                aps[i] = aps[i] - a * orb.bg[i][j];
                aqs[i] = aqs[i] - a * orb.bp[i][j];
            }
        }
    }
}

} // namespace feff::fovrg
