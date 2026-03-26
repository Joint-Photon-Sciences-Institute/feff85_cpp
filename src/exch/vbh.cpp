#include "vbh.hpp"
#include <cmath>

namespace feff::exch {

// gamma = 4/3 * asm / (1 - asm) where asm = 2^(-1/3)
// Precomputed to full precision as in Fortran
static constexpr double gamma_vbh = 5.129762802484098;

double flarge(double x) {
    return (1.0 + x * x * x) * std::log(1.0 + 1.0 / x) + x / 2.0 - x * x - 1.0 / 3.0;
}

void vbh(double rs, double xmag, double& vxc) {
    vxc = 0.0;
    if (rs > 1000.0) {
        // Transform to code units (Hartrees) from Rydbergs
        vxc = vxc / 2.0;
        return;
    }

    double epc  = -0.0504 * flarge(rs / 30.0);
    double efc  = -0.0254 * flarge(rs / 75.0);
    double xmup = -0.0504 * std::log(1.0 + 30.0 / rs);

    double vu = gamma_vbh * (efc - epc);

    double alg = -1.22177412 / rs + vu;
    double blg = xmup - vu;
    vxc = alg * std::pow(xmag, 1.0 / 3.0) + blg;

    // Transform to code units (Hartrees) from Rydbergs
    vxc = vxc / 2.0;
}

} // namespace feff::exch
