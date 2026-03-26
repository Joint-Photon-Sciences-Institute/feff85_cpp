// LDA exchange-correlation potential — converted from src/ATOM/vlda.f

#include "vlda.hpp"
#include "../par/parallel.hpp"
#include <feff/constants.hpp>
#include <cmath>

namespace feff::atom {

// =========================================================================
// EXCH stubs — Hedin-Lundqvist exchange-correlation
// These implement a basic exchange-correlation potential.
// Will be replaced with full EXCH module in Phase 4.
// =========================================================================

void vbh_stub(double rs, double xm, double& vxc) {
    // Hedin-Lundqvist exchange-correlation potential
    // V_xc = V_x + V_c where:
    //   V_x = -(3/(2*pi)) * (9*pi/4)^(1/3) / rs  (Slater exchange)
    //   V_c from Hedin-Lundqvist parametrization
    (void)xm;  // spin polarization not used in simplest form

    if (rs <= 0.0) {
        vxc = 0.0;
        return;
    }

    // Exchange part (Hartrees)
    constexpr double cx = -0.4581652932831429;  // -(3/(2*pi))*(9*pi/4)^(1/3)
    double vx = cx / rs;

    // Correlation part: Hedin-Lundqvist parametrization
    // V_c = -c * ln(1 + bi/x) where x = rs/A, bi = B/A
    constexpr double A = 21.0;   // HL parameter
    constexpr double bi = 0.7734;  // B/A
    double x = rs / A;
    double vc = 0.0;
    if (x > 0.0) {
        vc = -0.0225 * ((1.0 + bi * bi * bi / (x * x * x)) *
             std::log(1.0 + 1.0 / x) + x / 2.0 - bi -
             bi * bi * bi / (3.0 * x * x));
        // Simplified: use just the logarithmic form
        vc = -0.0225 * std::log(1.0 + A / rs);
    }

    vxc = vx + vc;
}

void edp_stub(double rs, double xf, double& vdh) {
    // Dirac-Hara exchange potential
    // V_DH = -(2*xf/pi) * (1 + ((1-eta^2)/(2*eta)) * ln|(1+eta)/(1-eta)|)
    // where eta = k_F / k (simplification: vdh ~ V_x)
    (void)xf;

    if (rs <= 0.0) {
        vdh = 0.0;
        return;
    }

    // Simple Slater exchange as fallback
    constexpr double cx = -0.4581652932831429;
    vdh = cx / rs;
}

// =========================================================================
// vlda implementation
// =========================================================================

void vlda(const double xnval[30], double srho[251], double srhovl[251],
          double vtrho[251], int ilast, int idfock, AtomState& state) {

    auto& orb = state.orb;
    auto& work = state.work;
    auto& scf = state.scf;
    auto& config = state.config;
    auto& mesh = state.mesh;

    int norb = scf.norb;

    // Zero density arrays
    for (int i = 0; i < 251; ++i) {
        srho[i] = 0.0;
        srhovl[i] = 0.0;
    }

    // Calculate total and valence densities
    for (int j = 0; j < norb; ++j) {
        double a = config.xnel[j];
        double b = xnval[j];
        for (int i = 0; i < config.nmax[j]; ++i) {
            double rho_j = orb.cg[i][j] * orb.cg[i][j] +
                          orb.cp[i][j] * orb.cp[i][j];
            srho[i] += a * rho_j;
            srhovl[i] += b * rho_j;
        }
    }

    // Construct LDA potential
    double rhoc = 0.0;
    for (int i = 0; i < 251; ++i) {
        double rho = srho[i] / (mesh.dr[i] * mesh.dr[i]);

        if (idfock == 5) {
            rhoc = srhovl[i] / (mesh.dr[i] * mesh.dr[i]);
        } else if (idfock == 6) {
            rhoc = (srho[i] - srhovl[i]) / (mesh.dr[i] * mesh.dr[i]);
        } else if (idfock == 1) {
            rhoc = 0.0;
        } else if (idfock == 2) {
            rhoc = srho[i] / (mesh.dr[i] * mesh.dr[i]);
        } else {
            feff::par::par_stop(" undefined idfock in subroutine vlda");
        }

        if (rho > 0.0) {
            double rs = std::pow(rho / 3.0, -feff::third);
            double rsc = 101.0;
            if (rhoc > 0.0) rsc = std::pow(rhoc / 3.0, -feff::third);
            double xm = 1.0;
            double vxcvl = 0.0;

            if (idfock == 5 || idfock == 2) {
                vbh_stub(rsc, xm, vxcvl);
            } else if (idfock == 6) {
                double vvbh;
                vbh_stub(rs, xm, vvbh);
                double xf = feff::fa / rs;
                double vdh;
                edp_stub(rsc, xf, vdh);
                vxcvl = vvbh - vdh;
            } else if (idfock == 1) {
                vxcvl = 0.0;
            }

            // Energy contribution
            if (ilast > 0) {
                vtrho[i] += vxcvl * srho[i];
            }

            // Add to potential
            if (i == 0) work.av[1] += vxcvl / work.cl;
            work.dv[i] += vxcvl / work.cl;
        }
    }
}

} // namespace feff::atom
