// Initialize orbital configuration and mesh for FOVRG.
// Converted from: src/FOVRG/inmuac.f

#include "inmuac.hpp"
#include "../../src/common/orbital_data.hpp"
#include <feff/dimensions.hpp>
#include <cmath>

namespace feff::fovrg {

void inmuac(int ihole, double xionin, int iunf, int ikap, FovrgState& state)
{
    auto& orb = state.orb;
    auto& config = state.config;
    auto& scf = state.scf;
    auto& lagrange = state.lagrange;
    auto& nuclear = state.nuclear;
    auto& mesh = state.mesh;
    int idim = mesh.idim;
    int nz = scf.nz;

    constexpr int nucm = 11;

    scf.testy = 1.0e-5;

    // Get orbital data for this element
    double xnval[30] = {};
    int iorb_dummy = 0;
    int iholep = 0;
    // getorb uses xmag as last arg; use en as dummy (same as Fortran)
    feff::common::getorb(nz, ihole, xionin, iunf,
                         scf.norb, scf.norbsc, iorb_dummy, iholep,
                         config.nq, config.kap, config.xnel, xnval, config.en);

    int norb = scf.norb;

    lagrange.ipl = 0;
    for (int i = 0; i < norb; i++) {
        config.en[i] = 0.0;
        lagrange.nre[i] = -1;
        int llq = std::abs(config.kap[i]);
        int l = llq + llq;

        // Find last tabulation point
        config.nmax[i] = 0;
        for (int j = idim - 1; j >= 0; j--) {
            if (std::abs(orb.cg[j][i]) >= 1.0e-11 ||
                std::abs(orb.cp[j][i]) >= 1.0e-11) {
                config.nmax[i] = j + 1;  // 1-based (matches Fortran convention)
                break;
            }
        }

        config.scc[i] = 0.3;
        if (config.xnel[i] < l) lagrange.nre[i] = 1;
        if (ikap == config.kap[i]) lagrange.ipl = lagrange.ipl + 1;
    }

    // Set up photoelectron orbital
    scf.norbsc = norb;
    scf.norb = norb + 1;
    config.xnel[norb] = 1.0;
    config.kap[norb] = ikap;
    config.nq[norb] = 9;

    // Nuclear radius
    nuclear.nuc = nucm;
}

} // namespace feff::fovrg
