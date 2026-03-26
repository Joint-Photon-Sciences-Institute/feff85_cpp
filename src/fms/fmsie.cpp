// Single energy FMS calculation.
// Converted from: fmsie.f
// Written by A. Ankudinov 06.1997, modified 2001-2002

#include "fmsie.hpp"
#include "fms_core.hpp"
#include "yprep.hpp"

#include <cmath>
#include <complex>
#include <vector>
#include <iostream>
#include <sstream>

namespace feff::fms {

void fmsie(bool verbose, int iph0, int nph, const int* lipotx,
           int ie, FeffComplex em, FeffComplex eref,
           const FeffComplex* ph,
           float rfms, int lfms, int nat,
           const int* iphat, const double* rath,
           Complexf* gtr,
           FMSData& data) {

    if (rfms <= 0.0f) return;

    // Default settings
    int minv = 0;        // LU decomposition
    float rdirec = 2.0f * rfms;
    float toler1 = 0.0f;
    float toler2 = 0.0f;

    // Convert double-precision coordinates to single
    std::vector<float> rat(nat * 3);
    for (int iat = 0; iat < nat; ++iat) {
        for (int j = 0; j < 3; ++j) {
            rat[iat * 3 + j] = static_cast<float>(rath[iat * 3 + j]);
        }
    }

    // Call yprep on first energy point or if lfms requests it
    int inclus = 0;
    static int prev_inclus = 0;
    if (ie == 1 || lfms == 0 || lfms == 2) {
        yprep(iph0, nat, inclus, iphat, rfms, rat.data(), data);
        prev_inclus = inclus;
    } else {
        inclus = prev_inclus;
    }

    if (inclus <= 1) return;

    // Print message on first energy
    if (ie == 1 && verbose) {
        std::cout << "        Doing FMS for a cluster of " << inclus
                  << " atoms around iph = " << iph0 << std::endl;
    }

    // Compute complex momentum
    FeffComplex dck = std::sqrt(2.0 * (em - eref));
    Complexf ck_val(static_cast<float>(dck.real()),
                    static_cast<float>(dck.imag()));
    Complexf ck_arr[nspx];
    ck_arr[0] = ck_val;

    // Convert phase shifts to single precision
    // Fortran: ph(1+abs(ill), ipp) with ill in [-lipotx, lipotx]
    // C++: ph[(lx+1) * iph + l] where l = abs(ill)
    // xphase layout: [isp][l+lx][iph] with l in [-lx, lx]
    std::vector<Complexf> xphase(nspx * (2 * lx + 1) * (nphx + 1), Complexf(0, 0));
    auto xph_idx = [](int isp, int l, int iph) -> int {
        return isp + nspx * ((l + lx) + (2 * lx + 1) * iph);
    };

    for (int ipp = 0; ipp <= nph; ++ipp) {
        for (int ill = -lipotx[ipp]; ill <= lipotx[ipp]; ++ill) {
            int l_abs = std::abs(ill);
            // ph is stored as ph[(lx+1) * iph + l]
            FeffComplex ph_val = ph[l_abs + (lx + 1) * ipp];
            xphase[xph_idx(0, ill, ipp)] = Complexf(
                static_cast<float>(ph_val.real()),
                static_cast<float>(ph_val.imag()));
        }
    }

    int iverb = 0;
    if (ie == 1) iverb = 1;
    if (!verbose) iverb = 0;

    int nsp = 1;
    int ispin = 0;
    bool lcalc[lx + 1];
    for (int i = 0; i <= lx; ++i) lcalc[i] = true;

    // Call FMS
    std::vector<Eigen::MatrixXcf> gg;
    fms(lfms, nsp, ispin, inclus, nph, ck_arr, lipotx, xphase.data(),
        ie, iverb, minv, rdirec, toler1, toler2, lcalc, gg, data);

    // Compute Green's function trace: gtr(il, ip) += Tr(gg) * exp(2i*delta)/(2l+1)
    Complexf coni_f(0.0f, 1.0f);
    constexpr int nsp_lx2 = nspx * (lx + 1) * (lx + 1);

    for (int ip = 0; ip <= nph; ++ip) {
        if (lfms == 0 && ip != iph0) continue;

        for (int il = 0; il <= lipotx[ip]; ++il) {
            int ix = il * il;  // Start index for this l block
            Complexf trace(0.0f, 0.0f);
            for (int im = 0; im < 2 * il + 1; ++im) {
                trace += gg[ip](ix + im, ix + im);
            }

            // gtr(il, ip) += trace * exp(2i * delta_l) / (2l+1)
            int gtr_idx = il + (lx + 1) * ip;
            gtr[gtr_idx] += trace *
                std::exp(2.0f * coni_f * xphase[xph_idx(0, il, ip)]) /
                static_cast<float>(2 * il + 1);
        }
    }
}

} // namespace feff::fms
