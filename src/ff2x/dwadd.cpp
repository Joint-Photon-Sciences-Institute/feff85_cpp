// Debye-Waller factor addition to paths.
// Converted from: src/FF2X/ff2gen.f (dwadd subroutine)
// Adds DW factors, cumulants, sums path contributions into cchi.

#include "dwadd.hpp"

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>

#include "../common/logging.hpp"

#include "../math/interpolation.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>

// Include actual Debye-Waller header
#include "../debye/debye.hpp"

namespace feff::ff2x {

static void pijump(double& phase, double phase_old) {
    while (phase - phase_old > pi) phase -= 2.0 * pi;
    while (phase - phase_old < -pi) phase += 2.0 * pi;
}

void dwadd(int ntotal, int nptot, const FF2xParams& p,
           const std::vector<PathListEntry>& path_list,
           FeffPadData& pad, const XsectData& xs,
           int nkx, const double xk0[], const double xkp[],
           FeffComplex cchi[], int iabs, int& nused) {

    auto& log = common::logger();
    nused = 0;

    bool dwcorr = (p.tk > 1.0e-3);
    int ne1 = xs.ne1;

    // Cumulant output file
    std::ofstream cum_out;
    if (p.alphat > 0.0) {
        cum_out.open("cum.dat");
        if (cum_out.is_open()) {
            cum_out << "# first and third cumulant for single scattering paths\n";
            cum_out << "# Einstein-Temp. =" << std::fixed << std::setprecision(2)
                    << p.thetae << "   alpha=" << std::setprecision(5) << p.alphat << "\n";
            cum_out << "#       file   sig1    sig2    sig3\n";
        }
    }

    // Cycle over all paths in the list
    for (int ilist = 0; ilist < ntotal; ++ilist) {
        // Find path index in feff.pad
        int ipath = -1;
        for (int j = 0; j < nptot; ++j) {
            if (path_list[ilist].ip == pad.index[j]) {
                ipath = j;
                break;
            }
        }
        if (ipath < 0) {
            char buf[128];
            std::snprintf(buf, sizeof(buf), " did not find path i=%d, ip=%d",
                          ilist, path_list[ilist].ip);
            log.wlog(buf);
            continue;
        }

        // Filter by curved-wave amplitude ratio
        if (pad.crit[ipath] < p.critcw) continue;

        // Compute Debye-Waller factor
        double sig2 = p.sig2g + path_list[ilist].sig2u;

        if (dwcorr && p.idwopt >= 0) {
            // Build temporary arrays for DW calculation (in Angstroms)
            int nl = pad.nleg[ipath];
            double rattmp[legtot + 1][3]; // (3, 0:legtot)
            int iztmp[legtot + 1];         // (0:legtot)

            for (int ileg = 0; ileg < nl; ++ileg) {
                iztmp[ileg + 1] = pad.iz[pad.ipot[ipath][ileg]];
                for (int j = 0; j < 3; ++j) {
                    rattmp[ileg + 1][j] = pad.rat[ipath][ileg][j] * bohr;
                }
            }
            // Central atom (index 0) = last leg
            iztmp[0] = iztmp[nl];
            for (int j = 0; j < 3; ++j) {
                rattmp[0][j] = rattmp[nl][j];
            }

            double sig2d = 0.0;
            double rs = pad.rnrmav;

            if (p.idwopt == 0) {
                // Correlated Debye model
                debye::sigms(p.tk, p.thetad, rs, legtot, nl,
                             rattmp, iztmp, sig2d);
            } else if (p.idwopt == 3) {
                // Classical Debye model
                debye::sigcl(p.tk, p.thetad, rs, legtot, nl,
                             rattmp, iztmp, sig2d);
            }
            // idwopt == 1 (EM) and 2 (RM) would call sigem/sigrm here

            sig2 += sig2d;
        }

        // Convert to code units
        sig2 /= (bohr * bohr);

        // First and third cumulants
        double sig1 = 0.0, sig3 = 0.0;
        if (p.alphat > zero && pad.nleg[ipath] == 2) {
            int nl = pad.nleg[ipath];
            int iz1 = pad.iz[pad.ipot[ipath][nl - 1]];
            int iz2 = pad.iz[pad.ipot[ipath][0]];

            if (p.thetae <= 0.0) {
                debye::sigte3(iz1, iz2, sig2, p.alphat, p.thetad,
                              pad.reff[ipath], sig1, sig3);
            } else {
                debye::sigm3(sig1, sig2, sig3, p.tk, p.alphat, p.thetae);
            }

            if (cum_out.is_open()) {
                char buf[64];
                std::snprintf(buf, sizeof(buf), "%10d%9.5f%9.5f %9.7f",
                              pad.index[ipath], sig1 * bohr,
                              sig2 * bohr * bohr, sig3 * bohr * bohr * bohr);
                cum_out << buf << "\n";
            }
        }

        // Apply DW factor and cumulants to achi/phchi
        // Use effective s02 from xsect (handles S02=0 meaning "use calculated")
        double s02_use = xs.s02_eff;
        if (p.mbconv > 0) s02_use = 1.0;

        for (int i = 0; i < ne1; ++i) {
            FeffComplex ck_d(pad.ck[i].real(), pad.ck[i].imag());
            FeffComplex dw = std::exp(-2.0 * sig2 * ck_d * ck_d);
            FeffComplex dw1 = std::exp(2.0 * coni * ck_d * sig1);
            FeffComplex dw3 = std::exp((-4.0 * coni * ck_d * ck_d * ck_d * sig3) / 3.0);
            dw *= dw1 * dw3;

            double phdw = 0.0;
            if (std::abs(dw) > 0.0)
                phdw = std::atan2(dw.imag(), dw.real());

            pad.achi[ipath][i] *= static_cast<float>(
                std::abs(dw) * s02_use * pad.deg[ipath]);
            pad.phchi[ipath][i] += static_cast<float>(phdw);
        }

        // Remove 2*pi jumps in phase
        for (int i = 1; i < ne1; ++i) {
            double curr = pad.phchi[ipath][i];
            double old_val = pad.phchi[ipath][i - 1];
            pijump(curr, old_val);
            pad.phchi[ipath][i] = static_cast<float>(curr);
        }

        // Interpolate and sum into cchi
        FeffComplex ccpath[nfinex];
        for (int ik = 0; ik < nkx; ++ik) {
            double achi0 = 0.0, phchi0 = 0.0;
            // Linear interpolation (terp1) on single-precision arrays
            math::terp1(pad.xk.data(), pad.achi[ipath].data(), ne1, xk0[ik], achi0);
            math::terp1(pad.xk.data(), pad.phchi[ipath].data(), ne1, xk0[ik], phchi0);

            ccpath[ik] = achi0 * std::exp(
                coni * (2.0 * xk0[ik] * static_cast<double>(pad.reff[ipath]) + phchi0));
            cchi[ik] += ccpath[ik];
        }
        nused++;
    }

    if (cum_out.is_open()) cum_out.close();
}

} // namespace feff::ff2x
