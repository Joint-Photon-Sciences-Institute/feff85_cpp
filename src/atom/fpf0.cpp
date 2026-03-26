// Atomic scattering amplitude f0 and oscillator strengths.
// Converted from src/ATOM/fpf0.f
//
// Writes fpf0.dat containing:
//   - Total energy contribution to f'
//   - Oscillator strengths for core-to-core dipole transitions
//   - f0(Q) scattering amplitude on a grid of 0.5 Angstrom^-1 steps
//
// JSON output from the Fortran version is omitted (json_module dependency).

#include "fpf0.hpp"
#include "../math/sommerfeld.hpp"
#include "../par/parallel.hpp"
#include "../common/logging.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <cmath>
#include <cstdio>
#include <fstream>

namespace feff::atom {

void fpf0(int iz, int iholep, const double srho[251], const double dr[251],
           double hx, const double dgc0[251], const double dpc0[251],
           const double dgc[], const double dpc[],
           double eatom, const double xnel[30], int norb,
           const double eorb[30], const int kappa[30])
{
    // dgc/dpc are flat arrays with layout dgc(251, 30, 0:nphx)
    // Access: dgc[i + 251*(j + 30*iph)]  (all 0-based)
    // For iph=0: dgc[i + 251*j]
    constexpr int np = 251;
    constexpr int stride_j = 251;
    constexpr int stride_iph = 251 * 30;

    double xpc[251] = {};
    double xqc[251] = {};

    // Output arrays
    double enosc[13] = {};
    double oscstr[13] = {};
    int index[13] = {};

    std::ofstream fout;
    bool open_16 = false;

    if (feff::par::state().master) {
        fout.open("fpf0.dat");
        open_16 = fout.is_open();

        double fpcorr = -std::pow(static_cast<double>(iz) / 82.5, 2.37);
        double etot_fp = eatom * feff::alphfs * feff::alphfs * 5.0 / 3.0;

        if (open_16) {
            char buf[256];
            std::snprintf(buf, sizeof(buf), " atom Z = %d", iz);
            fout << buf << "\n";
            std::snprintf(buf, sizeof(buf), "%19.5e%19.5e total energy part of fprime - 5/3*E_tot/mc**2",
                          etot_fp, fpcorr);
            fout << buf << "\n";
        }
    }

    // Get oscillator strengths
    enosc[0] = eorb[iholep];
    index[0] = iholep;
    int kinit = kappa[iholep];
    oscstr[0] = 2 * std::abs(kinit);
    int nosc = 1;

    double xmult1 = 0.0;
    double xmult2 = 0.0;

    for (int iorb = 0; iorb < norb; ++iorb) {
        if (xnel[iorb] > 0.0) {
            // Core orbital — check dipole selection rule
            int jkap = kappa[iorb];
            if (jkap + kinit == 0 || std::abs(jkap - kinit) == 1) {
                // Calculate reduced dipole matrix element
                int kdif = jkap - kinit;
                if (std::abs(kdif) > 1) kdif = 0;

                double twoj = 2.0 * std::abs(kinit) - 1.0;

                if (kdif == -1 && kinit > 0) {
                    xmult1 = 0.0;
                    xmult2 = std::sqrt(2.0 * (twoj + 1) * (twoj - 1) / twoj);
                } else if (kdif == -1 && kinit < 0) {
                    xmult1 = 0.0;
                    xmult2 = -std::sqrt(2.0 * (twoj + 1) * (twoj + 3) / (twoj + 2));
                } else if (kdif == 0 && kinit > 0) {
                    xmult1 = -std::sqrt((twoj + 1) * twoj / (twoj + 2));
                    xmult2 = -std::sqrt((twoj + 1) * (twoj + 2) / twoj);
                } else if (kdif == 0 && kinit < 0) {
                    xmult1 = std::sqrt((twoj + 1) * (twoj + 2) / twoj);
                    xmult2 = std::sqrt((twoj + 1) * twoj / (twoj + 2));
                } else if (kdif == 1 && kinit > 0) {
                    xmult1 = std::sqrt(2.0 * (twoj + 1) * (twoj + 3) / (twoj + 2));
                    xmult2 = 0.0;
                } else if (kdif == 1 && kinit < 0) {
                    xmult1 = -std::sqrt(2.0 * (twoj + 1) * (twoj - 1) / twoj);
                    xmult2 = 0.0;
                }

                double xk0 = std::abs(eorb[iorb] - eorb[iholep]) * feff::alphfs;

                for (int i = 0; i < np; ++i) {
                    double xj0 = std::sin(xk0 * dr[i]) / (xk0 * dr[i]);
                    // dgc(i, iorb, 0) = dgc[i + 251*iorb]
                    // dpc(i, iorb, 0) = dpc[i + 251*iorb]
                    xpc[i] = (xmult1 * dgc0[i] * dpc[i + stride_j * iorb] +
                              xmult2 * dpc0[i] * dgc[i + stride_j * iorb]) * xj0;
                    xqc[i] = 0.0;
                }

                double xirf = 2.0;
                feff::math::somm(dr, xpc, xqc, hx, xirf, 0, np);

                oscstr[nosc] = xirf * xirf / 3.0;
                enosc[nosc] = eorb[iorb];
                index[nosc] = iorb;
                nosc++;
            }
        }
    }

    // Write oscillator strength information
    if (open_16) {
        char buf[256];
        std::snprintf(buf, sizeof(buf), " %d", nosc);
        fout << buf << "\n";
        for (int i = 0; i < nosc; ++i) {
            std::snprintf(buf, sizeof(buf), " %9.5f %12.3f %4d",
                          oscstr[i], enosc[i], index[i]);
            fout << buf << "\n";
        }
    }

    // Calculate and write f0(Q) on grid delq = 0.5 Angstrom^-1
    double dq = 0.5 * feff::bohr;  // convert to atomic units
    for (int iq = 0; iq < 81; ++iq) {
        double xk0 = dq * iq;
        // srho is 4*pi*density
        for (int i = 0; i < np; ++i) {
            double xj0 = 1.0;
            if (iq > 0) xj0 = std::sin(xk0 * dr[i]) / (xk0 * dr[i]);
            xpc[i] = srho[i] * (dr[i] * dr[i]) * xj0;
            xqc[i] = 0.0;
        }
        double xirf = 2.0;
        feff::math::somm(dr, xpc, xqc, hx, xirf, 0, np);

        if (open_16) {
            char buf[64];
            std::snprintf(buf, sizeof(buf), " %5.1f %9.4f",
                          0.5 * iq, xirf);
            fout << buf << "\n";
        }
    }

    if (open_16) {
        fout.close();
    }
}

} // namespace feff::atom
