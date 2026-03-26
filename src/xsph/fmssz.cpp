// FMS self-consistency integration for szlz calculations.
// Converted from src/XSPH/fmssz.f

#include "fmssz.hpp"
#include "acoef.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include "../common/logging.hpp"
#include <cmath>
#include <complex>
#include <cstring>

// Forward declarations for FMS routines
namespace feff::fms {
    void yprep(int iph0, int nat, int& inclus, const int iphat[],
               float rfms, const float rat[]);
    void fmsie(int lfms, int nsp, int ispin, int inclus, int nph,
               const std::complex<float> ck[], const int lipotx[],
               const std::complex<float>* xphase, int ie,
               int iverb, int minv, float rdirec, float toler1, float toler2,
               const bool lcalc[], std::complex<float>* gg);
}

namespace feff::xsph {

void fmssz(bool verbose, int iph0, int ie, FeffComplex em, FeffComplex eref,
           const FeffComplex* ph, int nph,
           float rfms, int lfms, int nat, const int iphat[],
           const double* rat, const float amat[], const int lipotx[],
           float* gctr, FeffComplex* gtr) {

    // ph indexed as ph[(lx+1) * iph + il]
    // gtr indexed as gtr[i1 + 2*(i2 + 2*(iop + 3*(il + (lx+1)*ip)))]
    constexpr int gtr_stride = 2 * 2 * 3 * (lx + 1);

    if (rfms <= 0.0f) {
        // No FMS, just accumulate central atom part
        for (int ip = 0; ip <= nph; ip++) {
            if (lfms != 0 || ip == iph0) {
                for (int lpp = 0; lpp <= lipotx[ip]; lpp++) {
                    int ix1 = lpp * lpp;
                    for (int im = 0; im < 2 * lpp + 1; im++) {
                        for (int iop = 0; iop < 3; iop++) {
                            for (int i2 = 0; i2 < 2; i2++) {
                                for (int i1 = 0; i1 < 2; i1++) {
                                    int g_idx = i1 + 2 * (i2 + 2 * (iop + 3 * (lpp + (lx + 1) * ip)));
                                    int a_idx = amat_index(im - lpp, i1, i2, iop, lpp);
                                    gctr[g_idx] += amat[a_idx];
                                }
                            }
                        }
                    }
                }
            }
        }
        return;
    }

    // Convert coordinates to single precision
    float ratf[3 * natx];
    for (int iat = 0; iat < nat; iat++) {
        for (int j = 0; j < 3; j++) {
            ratf[j + 3 * iat] = static_cast<float>(rat[j + 3 * iat]);
        }
    }

    int inclus = 0;
    if (ie == 0 || lfms == 0) {
        feff::fms::yprep(iph0, nat, inclus, iphat, rfms, ratf);
    }

    if (inclus <= 1) return;

    if (ie == 0 && verbose) {
        char slog[512];
        std::snprintf(slog, sizeof(slog),
                     "        Doing FMS for a cluster of %3d atoms around iph = %2d",
                     inclus, iph0);
        feff::common::logger().wlog(slog);
    }

    // Complex momentum
    FeffComplex dck = std::sqrt(2.0 * (em - eref));
    std::complex<float> ck_f(static_cast<float>(std::real(dck)),
                              static_cast<float>(std::imag(dck)));

    // Phase shifts (single precision)
    std::complex<float> xphase[nspx * (2 * lx + 1) * (nphx + 1)];
    for (int ipp = 0; ipp <= nph; ipp++) {
        for (int ill = -lipotx[ipp]; ill <= lipotx[ipp]; ill++) {
            int ph_idx = std::abs(ill) + (lx + 1) * ipp;
            float rp = static_cast<float>(std::real(ph[ph_idx]));
            float ip_val = static_cast<float>(std::imag(ph[ph_idx]));
            // Index into xphase: xphase[0][ill+lx][ipp]
            xphase[(ill + lx) + (2 * lx + 1) * ipp] = std::complex<float>(rp, ip_val);
        }
    }

    int iverb = (ie == 0) ? 1 : 0;
    int nsp = 1, ispin_fms = 0, minv = 0;
    float rdirec = 2.0f * rfms, toler1 = 0.0f, toler2 = 0.0f;
    bool lcalc[lx + 1];
    for (int i = 0; i <= lx; i++) lcalc[i] = true;

    constexpr int gg_dim = nspx * (lx + 1) * (lx + 1);
    std::complex<float> gg[gg_dim * gg_dim * (nphx + 1)];
    std::memset(gg, 0, sizeof(gg));

    feff::fms::fmsie(lfms, nsp, ispin_fms, inclus, nph, &ck_f, lipotx, xphase,
                     ie, iverb, minv, rdirec, toler1, toler2, lcalc, gg);

    // Accumulate results
    for (int ip = 0; ip <= nph; ip++) {
        if (lfms != 0 || ip == iph0) {
            for (int lpp = 0; lpp <= lipotx[ip]; lpp++) {
                int ix1 = lpp * lpp;
                for (int im = 0; im < 2 * lpp + 1; im++) {
                    for (int iop = 0; iop < 3; iop++) {
                        for (int i2 = 0; i2 < 2; i2++) {
                            for (int i1 = 0; i1 < 2; i1++) {
                                int a_idx = amat_index(im - lpp, i1, i2, iop, lpp);
                                int g_idx = i1 + 2 * (i2 + 2 * (iop + 3 * (lpp + (lx + 1) * ip)));
                                int gg_idx = (ix1 + im) + gg_dim * ((ix1 + im) + gg_dim * ip);
                                auto gg_val = gg[gg_idx];

                                gtr[g_idx] += FeffComplex(amat[a_idx], 0.0) *
                                    FeffComplex(gg_val.real(), gg_val.imag());
                                gctr[g_idx] += amat[a_idx];
                            }
                        }
                    }
                }
            }
        }
    }
}

} // namespace feff::xsph
