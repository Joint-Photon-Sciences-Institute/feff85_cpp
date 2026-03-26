// Initialize orbital occupations and SCF parameters.
// Converted from: src/ATOM/inmuat.f
//
// INDEX CONVENTION: Loop indices remain 1-based as in Fortran.
// Array accesses use F() macro for 0-based C++ arrays.

#include "inmuat.hpp"
#include "../common/orbital_data.hpp"
#include "../par/parallel.hpp"
#include <cmath>

namespace feff::atom {

#define F(arr, i) (arr)[(i) - 1]

void inmuat(int ihole, double xionin, int iunf, double xnval[],
            int& iholep, double xmag[], int& iorb, AtomState& state)
{
    auto& config  = state.config;
    auto& scf     = state.scf;
    auto& lagrange = state.lagrange;
    auto& nuclear = state.nuclear;
    auto& mesh    = state.mesh;

    // Constants from Fortran DATA statement
    constexpr int nucm = 11;
    constexpr int nesn = 50;
    constexpr int ideps = 435;

    mesh.ndor = 10;
    scf.testy = 1.0e-05;
    scf.teste = 5.0e-06;
    scf.rap[0] = 100.0;  // rap(1)
    scf.rap[1] = 10.0;   // rap(2)

    // Zero out arrays
    for (int i = 1; i <= 30; ++i) {
        F(config.en, i) = 0.0;
        F(xmag, i) = 0.0;
        F(xnval, i) = 0.0;
    }

    // Get orbital data
    feff::common::getorb(scf.nz, ihole, xionin, iunf,
                         scf.norb, scf.norbsc, iorb, iholep,
                         config.nq, config.kap, config.xnel,
                         xnval, xmag);

    // Check total electron count
    double xk = 0.0;
    for (int i = 1; i <= scf.norb; ++i) {
        xk += F(config.xnel, i);
    }
    if (std::abs(scf.nz - xionin - xk) > 0.001) {
        feff::par::par_stop("check number of electrons in getorb.f");
    }
    scf.norbsc = scf.norb;

    mesh.nes = nesn;
    nuclear.nuc = nucm;

    // Zero Lagrange parameters
    for (int i = 1; i <= ideps; ++i) {
        F(lagrange.eps, i) = 0.0;
    }

    mesh.idim = 251;
    if (mesh.idim % 2 == 0) mesh.idim = mesh.idim - 1;

    lagrange.ipl = 0;

    for (int i = 1; i <= scf.norb; ++i) {
        F(lagrange.nre, i) = -1;
        int llq = std::abs(F(config.kap, i));
        int l = llq + llq;
        if (F(config.kap, i) < 0) llq = llq - 1;
        if (llq < 0 || llq >= F(config.nq, i) || llq > 3) {
            feff::par::par_stop("kappa out of range, check getorb.f");
        }
        F(config.nmax, i) = mesh.idim;
        F(config.scc, i) = 0.3;
        if (F(config.xnel, i) < l) F(lagrange.nre, i) = 1;
        if (F(config.xnel, i) < 0.5) F(config.scc, i) = 1.0;
        for (int j = 1; j <= i - 1; ++j) {
            if (F(config.kap, j) == F(config.kap, i)) {
                if (F(lagrange.nre, j) > 0 || F(lagrange.nre, i) > 0)
                    lagrange.ipl = lagrange.ipl + 1;
            }
        }
    }
}

#undef F

} // namespace feff::atom
