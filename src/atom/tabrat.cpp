// Orbital statistics tabulation — converted from src/ATOM/tabrat.f
//
// Tabulates orbital results: number of points, average values of r^n,
// and overlap integrals between orbitals of the same symmetry.
// Output goes to the logger when print level is sufficient.

#include "tabrat.hpp"
#include "radial_integrals.hpp"
#include "../common/logging.hpp"
#include <feff/constants.hpp>
#include <cmath>
#include <cstdio>

namespace feff::atom {

void tabrat(AtomState& state) {
    auto& config = state.config;
    auto& scf    = state.scf;
    auto& orb    = state.orb;
    auto& work   = state.work;
    auto& mesh   = state.mesh;

    int norb = scf.norb;

    // Orbital identification strings
    // Fortran ttire: 's ', 'p*', 'p ', 'd*', 'd ', 'f*', 'f ', 'g*', 'g '
    const char* ttire[9] = {"s ", "p*", "p ", "d*", "d ", "f*", "f ", "g*", "g "};
    char titre[30][3] = {};

    for (int i = 0; i < norb; ++i) {
        int j;
        if (config.kap[i] > 0) {
            j = 2 * config.kap[i];
        } else {
            j = -2 * config.kap[i] - 1;
        }
        // j is 1-based index into ttire in Fortran; 0-based here
        std::snprintf(titre[i], 3, "%s", ttire[j - 1]);
    }

    // Compute mbi array for r^n powers: n = 6,4,2,1,-1,-2,-3
    // Fortran: do i=2,8; mbi(i) = 8-i - i/3 - i/4 + i/8
    int mbi[9] = {};
    for (int i = 2; i <= 8; ++i) {
        mbi[i] = 8 - i - i / 3 - i / 4 + i / 8;
    }

    // Log header
    feff::common::logger().wlog(
        "number of electrons nel and average values of r**n in a.u.");

    char buf[256];
    std::snprintf(buf, sizeof(buf), "     nel     -E      n=%2d%10s%2d%10s%2d%10s%2d%10s%2d%10s%2d%10s%2d",
                  mbi[2], "", mbi[3], "", mbi[4], "", mbi[5], "",
                  mbi[6], "", mbi[7], "", mbi[8]);
    feff::common::logger().wlog(buf);

    // Tabulate average values of r^n for each orbital
    for (int i = 0; i < norb; ++i) {
        int llq = std::abs(config.kap[i]) - 1;
        int j = 8;
        if (llq <= 0) j = 7;

        double at[9] = {};
        for (int k = 2; k <= j; ++k) {
            at[k] = dsordf(i, i, mbi[k], 1, feff::zero,
                           orb, work, config, mesh);
        }

        // Format output line
        if (j == 8) {
            std::snprintf(buf, sizeof(buf), "%1d%s%6.3f%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e",
                          config.nq[i], titre[i], config.xnel[i],
                          -config.en[i] * feff::hart,
                          at[2], at[3], at[4], at[5], at[6], at[7]);
        } else {
            std::snprintf(buf, sizeof(buf), "%1d%s%6.3f%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e",
                          config.nq[i], titre[i], config.xnel[i],
                          -config.en[i] * feff::hart,
                          at[2], at[3], at[4], at[5], at[6]);
        }
        feff::common::logger().wlog(buf);
    }

    // Overlap integrals
    if (norb <= 1) return;

    feff::common::logger().wlog("          overlap integrals");

    for (int i = 0; i < norb - 1; ++i) {
        for (int j = i + 1; j < norb; ++j) {
            if (config.kap[j] == config.kap[i]) {
                double ovlp = dsordf(i, j, 0, 1, feff::zero,
                                     orb, work, config, mesh);
                std::snprintf(buf, sizeof(buf), "    %3d%s%3d%s%14.7f",
                              config.nq[i], titre[i], config.nq[j], titre[j], ovlp);
                feff::common::logger().wlog(buf);
            }
        }
    }
}

} // namespace feff::atom
