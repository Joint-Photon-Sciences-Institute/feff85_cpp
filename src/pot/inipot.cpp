// Initialize potential arrays to zero.
// Converted from src/POT/inipot.f

#include "inipot.hpp"
#include <cstring>

namespace feff::pot {

void inipot(double* dgc, double* dpc,
            double* edenvl, double* vvalgs, double* xnmues)
{
    // dgc, dpc: dimension (251, 30, 0:nphx+1)
    // Total elements: 251 * 30 * (nphx+2)
    constexpr int dgc_size = 251 * 30 * (nphx + 2);
    std::memset(dgc, 0, dgc_size * sizeof(double));
    std::memset(dpc, 0, dgc_size * sizeof(double));

    // edenvl, vvalgs: dimension (251, 0:nphx)
    // Total elements: 251 * (nphx+1)
    constexpr int ev_size = 251 * (nphx + 1);
    std::memset(edenvl, 0, ev_size * sizeof(double));
    std::memset(vvalgs, 0, ev_size * sizeof(double));

    // xnmues: dimension (0:lx, 0:nphx)
    // Total elements: (lx+1) * (nphx+1)
    constexpr int xn_size = (lx + 1) * (nphx + 1);
    std::memset(xnmues, 0, xn_size * sizeof(double));
}

} // namespace feff::pot
