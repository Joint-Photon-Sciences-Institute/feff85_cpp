// Full FMS orchestrator.
// Converted from: fmstot.f
// Drives the energy loop for XANES/FMS calculations.

#include "fmstot.hpp"
#include "fms_core.hpp"
#include "xprep.hpp"

#include <cmath>
#include <complex>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>

#include "../common/file_io.hpp"
#include "../common/logging.hpp"
#include "../common/pad_io.hpp"
#include "../math/bcoef.hpp"
#include "../math/wigner.hpp"

namespace feff::fms {

void fmstot(float rclust, int idwopt, double tk, double thetad, double sigma2,
            int* lmaxph, int nat, const int* iphat, const double* ratdbl,
            int ipol, int ispin, int le2, double angks,
            const FeffComplex ptz[3][3],
            int minv, float rdirec, float toler1, float toler2) {

    constexpr int npadx = 8;

    // Read phase.pad
    auto phdata = feff::common::read_xsph("phase.pad");
    int ne   = phdata.ne;
    int ne1  = phdata.ne1;
    int ne3  = phdata.ne3;
    int nph  = phdata.nph;
    int ihole = phdata.ihole;

    // Green's function trace accumulated over energies
    std::vector<Complexf> gtr(ne, Complexf(0.0f, 0.0f));

    if (rclust > 0.0f) {
        std::cout << " Number of energy points = " << ne << std::endl;

        // Convert coordinates to single precision (Bohr)
        std::vector<float> rat(nat * 3);
        for (int iat = 0; iat < nat; ++iat) {
            for (int j = 0; j < 3; ++j) {
                rat[iat * 3 + j] = static_cast<float>(ratdbl[iat * 3 + j]);
            }
        }

        float rnrmax = static_cast<float>(phdata.rnrmav);
        float temper = static_cast<float>(tk);
        float thetax = static_cast<float>(thetad);
        float sig2   = static_cast<float>(sigma2);

        // Prepare cluster geometry
        FMSData data;
        int inclus = 0;
        int iph0 = 0;
        int npot = nph;

        xprep(iph0, idwopt, nat, inclus, npot, iphat, rclust,
              rat.data(), phdata.iz.data(), rnrmax, temper, thetax, sig2,
              minv, rdirec, data);

        if (inclus > 1) {
            std::cout << " Doing FMS for a cluster of " << inclus
                      << " atoms around iph = " << iph0 << std::endl;
            std::cout << " Please, wait (updates every 20 points) ..." << std::endl;

            int nsp = 1;
            if (std::abs(ispin) == 1) nsp = nspx;

            // Compute bcoef matrices
            // setkap equivalent: determine kinit from ihole
            int kinit = -(ihole + 1);  // simplified; full mapping in pipeline
            int linit = (std::abs(kinit) - 1 + (kinit > 0 ? 1 : 0));
            (void)linit;

            int kind[8], lind[8];
            std::vector<FeffComplex> bmat0((2 * lx + 1) * 2 * 8 * (2 * lx + 1) * 2 * 8);
            feff::math::bcoef(kinit, ipol, ptz, le2, false, ispin, angks,
                              kind, lind, bmat0.data());

            // Phase shifts and momentum at each energy
            std::vector<Complexf> xphase(nspx * (2 * lx + 1) * (nphx + 1));
            Complexf ck[nspx];

            // lcalc: which l-channels to compute
            bool lcalc[lx + 1];
            for (int i = 0; i <= lx; ++i) lcalc[i] = false;
            for (int k1 = 0; k1 < 8; ++k1) {
                if (lind[k1] >= 0 && lind[k1] <= lx)
                    lcalc[lind[k1]] = true;
            }

            // Energy loop
            for (int ie = 0; ie < ne; ++ie) {
                // Compute complex momentum for each spin
                for (int isp = 0; isp < nsp; ++isp) {
                    FeffComplex dck = std::sqrt(
                        2.0 * (phdata.em[ie] - phdata.eref[ie * nsp + isp]));
                    ck[isp] = Complexf(static_cast<float>(dck.real()),
                                       static_cast<float>(dck.imag()));
                }

                // Convert phase shifts to single precision
                std::fill(xphase.begin(), xphase.end(), Complexf(0, 0));
                for (int ipp = 0; ipp <= nph; ++ipp) {
                    for (int isp = 0; isp < nsp; ++isp) {
                        for (int ill = -lmaxph[ipp]; ill <= lmaxph[ipp]; ++ill) {
                            int ph_idx = ie + ne * ((ill + ltot) +
                                (2 * ltot + 1) * (isp + nspx * ipp));
                            FeffComplex ph_val = phdata.ph[ph_idx];
                            int xph_idx = isp + nspx * ((ill + lx) + (2 * lx + 1) * ipp);
                            xphase[xph_idx] = Complexf(
                                static_cast<float>(ph_val.real()),
                                static_cast<float>(ph_val.imag()));
                        }
                    }
                }

                int iverb = 0;
                if ((ie + 1) % 20 == 0) iverb = 1;

                // Call FMS
                std::vector<Eigen::MatrixXcf> gg;
                fms(0, nsp, ispin, inclus, npot, ck, lmaxph, xphase.data(),
                    ie + 1, iverb, minv, rdirec, toler1, toler2, lcalc,
                    gg, data);

                // Accumulate gtr from bmat and gg
                for (int k1 = 0; k1 < 8; ++k1) {
                    for (int is1 = 0; is1 < nsp; ++is1) {
                        for (int k2 = 0; k2 < 8; ++k2) {
                            for (int is2 = 0; is2 < nsp; ++is2) {
                                if (lind[k2] < 0 || lind[k1] < 0) continue;

                                int ix1 = nsp * (lind[k1] * lind[k1] + lind[k1]);
                                int ix2 = nsp * (lind[k2] * lind[k2] + lind[k2]);
                                int ms1 = is1;
                                int ms2 = is2;

                                for (int m1 = -lind[k1]; m1 <= lind[k1]; ++m1) {
                                    for (int m2 = -lind[k2]; m2 <= lind[k2]; ++m2) {
                                        int bm_idx = feff::math::bmat_index(
                                            m2, ms2, k2, m1, ms1, k1);

                                        int gg_r = ix1 + nsp * m1 + is1;
                                        int gg_c = ix2 + nsp * m2 + is2;

                                        FeffComplex rkk1 = phdata.rkk[ie + ne * (k1 + 8 * is1)];
                                        FeffComplex rkk2 = phdata.rkk[ie + ne * (k2 + 8 * is2)];

                                        Complexf gg_val = gg[0](gg_r, gg_c);
                                        Complexf bm_val(
                                            static_cast<float>(bmat0[bm_idx].real()),
                                            static_cast<float>(bmat0[bm_idx].imag()));

                                        gtr[ie] += gg_val * bm_val *
                                            Complexf(static_cast<float>((rkk1 * rkk2).real()),
                                                     static_cast<float>((rkk1 * rkk2).imag()));
                                    }
                                }
                            }
                        }
                    }
                }
            } // energy loop
        } // inclus > 1
    } // rclust > 0

    // Write fms.bin
    std::ofstream fout("fms.bin");
    if (fout.is_open()) {
        fout << "FMS rfms=" << rclust * static_cast<float>(feff::bohr) << "\n";
        fout << ne << " " << ne1 << " " << ne3 << " " << nph << " " << npadx << " 1\n";

        std::vector<FeffComplex> dum(ne);
        for (int ie = 0; ie < ne; ++ie) {
            dum[ie] = FeffComplex(
                static_cast<double>(gtr[ie].real()),
                static_cast<double>(gtr[ie].imag()));
        }
        feff::common::write_pad_complex(fout, npadx, dum.data(), ne);
        fout.close();
    }
}

} // namespace feff::fms
