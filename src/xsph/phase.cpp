// Phase shift calculation per unique potential.
// Converted from src/XSPH/phase.f

#include "phase.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include "../exch/xcpot.hpp"
#include "../fovrg/dfovrg.hpp"
#include "../math/bessel.hpp"
#include "../math/phase_amplitude.hpp"
#include "../common/logging.hpp"
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <fstream>

namespace feff::xsph {

// Helper to index into ph: ph[ie][ll+ltot] where ll ranges [-ltot,ltot]
static inline int ph_idx(int ie, int ll) {
    return ie + nex * (ll + ltot);
}

void phase(int iph, double dx, double x0, const double ri[],
           int ne, int ne1, int ne3, const FeffComplex em[],
           int ixc, int nsp, int lreal, double rmt, double rnrm,
           double xmu, int iPl,
           const double vtot[], const double vvalgs[],
           const double edens[], const double dmag[], const double edenvl[],
           const double dgcn[][30], const double dpcn[][30],
           const double adgc[][30], const double adpc[][30],
           FeffComplex eref[], FeffComplex* ph, int& lmax,
           int iz, int ihole, double xion, int iunf,
           const double xnval[], int ispin) {

    // Workspace for xcpot
    double vxcrmu[nrptx], vxcimu[nrptx], gsrel[nrptx];
    double vvxcrm[nrptx], vvxcim[nrptx];
    FeffComplex p_arr[nrptx], q_arr[nrptx];
    FeffComplex v[nrptx], vval[nrptx];

    // Plasmon pole data
    double WpCorr[MxPole], AmpFac[MxPole];
    for (int i = 0; i < MxPole; i++) WpCorr[i] = -1.0e30;

    // Set imt and jri on Loucks grid
    int imt = static_cast<int>((std::log(rmt) + x0) / dx) + 1;
    int jri = imt + 1;
    int jri1 = jri + 1;
    if (jri1 > nrptx) {
        throw std::runtime_error("jri > nrptx in phase");
    }

    // Zero phase shifts and find xkmax
    double xkmax = 0.0;
    int ne12 = ne - ne3;
    for (int ie = 0; ie < ne; ie++) {
        for (int il = -ltot; il <= ltot; il++) {
            ph[ph_idx(ie, il)] = 0.0;
        }
        if (ie < ne12 && xkmax < std::real(em[ie])) {
            xkmax = std::real(em[ie]);
        }
    }
    xkmax = std::sqrt(xkmax * 2.0);

    // Use kmax to find lmax
    double prefac = 0.7;
    lmax = static_cast<int>(prefac * rmt * xkmax);
    lmax = std::max(lmax, 5);
    if (lmax > ltot) {
        int ik = static_cast<int>(std::round(ltot / rmt / bohr / prefac));
        char slog[512];
        std::snprintf(slog, sizeof(slog),
                     "      Phase shift calculation is accurate to k=%d", ik);
        feff::common::logger().wlog(slog);
        feff::common::logger().wlog("      See FEFF document to increase the range.");
    }
    lmax = std::min(lmax, ltot);

    // Read plasmon pole data if needed
    if (iPl > 0 && ixc == 0) {
        std::ifstream fin("exc.dat");
        if (fin.is_open()) {
            std::string line;
            int ipole = 0;
            while (ipole < MxPole && std::getline(fin, line)) {
                if (line.empty() || line[0] == '#' || line[0] == '*' ||
                    line[0] == 'c' || line[0] == 'C') continue;
                double wp_val, gam_val, amp_val;
                if (std::sscanf(line.c_str(), "%lf %lf %lf", &wp_val, &gam_val, &amp_val) == 3) {
                    double rs_local = std::pow(3.0 / (4.0 * pi * edens[jri + 1]), third);
                    double wp_calc = std::sqrt(3.0 / (rs_local * rs_local * rs_local));
                    WpCorr[ipole] = (wp_val / hart) / wp_calc;
                    AmpFac[ipole] = amp_val;
                    ipole++;
                }
            }
            if (ipole < MxPole) WpCorr[ipole] = -1.0e30;
        }
    }

    int ifirst = 0;
    int index = ixc;
    feff::atom::FovrgState state;

    // Calculate phase shifts for each energy
    for (int ie = 0; ie < ne12; ie++) {
        feff::exch::xcpot(iph, ie, index, lreal, ifirst, jri,
                          em[ie], xmu,
                          vtot, vvalgs, edens, dmag, edenvl,
                          eref[ie], v, vval, iPl, WpCorr, AmpFac,
                          vxcrmu, vxcimu, gsrel, vvxcrm, vvxcim, rnrm);

        if (std::real(em[ie]) < -10.0 || std::real(em[ie]) > 3.0e2) continue;

        // Complex momentum squared referenced to energy-dep XC
        FeffComplex p2 = em[ie] - eref[ie];
        if (lreal > 1 && ie < ne1) p2 = std::real(p2);
        FeffComplex ck = std::sqrt(2.0 * p2 + (p2 * alphfs) * (p2 * alphfs));
        FeffComplex xkmt = rmt * ck;
        if (std::real(p2) <= 0.0 && std::imag(p2) <= 0.0) continue;

        // Compute Bessel functions at muffin-tin
        FeffComplex jl[ltot + 2], nl[ltot + 2];
        feff::math::besjn(xkmt, jl, nl);

        int ncycle = (ixc % 10 < 5) ? 0 : 3;

        for (int ll = -lmax; ll <= lmax; ll++) {
            int il = std::abs(ll) + 1;
            // Nonlocal exchange unstable for high il
            int nc = ncycle;
            if (il * dx > 0.50) nc = 0;

            int ikap = ll - 1;
            if (ll > 0) ikap = ll;
            int ilp = il + 1;
            if (ikap > 0) ilp = il - 1;
            int ic3 = 0;

            if (nsp == 1 && ispin == 0) {
                // Remove spin-orbit interaction
                if (ll != 0) ic3 = 1;
                ikap = -1 - std::abs(ll);
                ilp = il + 1;
            }

            // Solve Dirac equation
            int irr = -1;
            FeffComplex pu, qu;
            FeffComplex p2_local = p2;
            int jlast = jri;
            feff::fovrg::dfovrg(nc, ikap, rmt, jlast, jri, p2_local, dx,
                               ri, v, vval, dgcn, dpcn, adgc, adpc,
                               xnval, pu, qu, p_arr, q_arr,
                               iz, ihole, xion, iunf, irr, ic3, state);

            // Extract phase shift from boundary matching
            FeffComplex temp;
            // il-1 and ilp-1 for 0-based indexing of jl, nl arrays
            feff::math::phamp(rmt, pu, qu, ck, jl[il - 1], nl[il - 1],
                             jl[ilp - 1], nl[ilp - 1], ikap,
                             ph[ph_idx(ie, ll)], temp);

            // Cut if phase shifts become too small
            if (std::abs(ph[ph_idx(ie, ll)]) < 1.0e-6 && ll >= 4) break;
            // Rivas cut
            if (std::abs(std::exp(FeffComplex(0, 2) * ph[ph_idx(ie, ll)]) - 1.0) < 1.0e-5) {
                ph[ph_idx(ie, ll)] = 0.0;
            }
            if (std::abs(ph[ph_idx(ie, ll)]) < 1.0e-5 && ll >= 4) break;
        }
    }

    // Set eref for remaining points
    for (int ie = ne12; ie < ne; ie++) {
        eref[ie] = eref[ne1 - 1];
    }
}

} // namespace feff::xsph
