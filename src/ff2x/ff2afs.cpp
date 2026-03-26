// Anomalous scattering f' / DANES calculation.
// Converted from: src/FF2X/ff2afs.f
// Calculate anomalous scattering amplitude for a given edge.

#include "ff2afs.hpp"
#include "ff2x.hpp"
#include "rdfbin.hpp"
#include "rdxbin.hpp"
#include "feffdt.hpp"
#include "exconv.hpp"
#include "fprime.hpp"
#include "xscorr.hpp"
#include "dwadd.hpp"

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>

#include "../common/logging.hpp"
#include "../common/pad_io.hpp"

#include "../math/interpolation.hpp"
#include "../math/convolution.hpp"

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

namespace feff::ff2x {

void ff2afs(const FF2xParams& p, int iabs) {
    auto& log = common::logger();

    // Read FMS Green's function from fms.bin (same as ff2xmu)
    int nip = p.ipmax - p.ipmin + 1;
    std::vector<std::vector<FeffComplex>> gtrful(nip);
    int ne_fms = 0;
    int ntfms = 0;
    std::string titfms;

    {
        std::ifstream fms_in("fms.bin");
        if (fms_in.is_open()) {
            ntfms = 1;
            std::getline(fms_in, titfms);
            int ne1_f = 0, ne3_f = 0, nph_f = 0, npadx_f = 0;
            {
                std::string line;
                std::getline(fms_in, line);
                std::istringstream iss(line);
                iss >> ne_fms >> ne1_f >> ne3_f >> nph_f >> npadx_f;
            }
            int total_pts = ne_fms * nip;
            std::vector<FeffComplex> gtrtemp(total_pts);
            common::read_pad_complex(fms_in, npadx_f, gtrtemp.data(), total_pts);
            int idx = 0;
            for (int iip = 0; iip < nip; ++iip) {
                gtrful[iip].resize(ne_fms);
                for (int j = 0; j < ne_fms; ++j)
                    gtrful[iip][j] = gtrtemp[idx + j];
                idx += ne_fms;
            }
        }
    }

    // Loop over spectrum indices
    for (int iip = p.ipmin; iip <= p.ipmax; iip += p.ipstep) {
        bool cross = is_cross_spectrum(iip);
        auto files = make_spectrum_files(iip);

        int iip_idx = iip - p.ipmin;
        std::vector<FeffComplex> gtr;
        if (ntfms > 0 && iip_idx < nip)
            gtr = gtrful[iip_idx];
        else
            gtr.resize(nex, FeffComplex(0.0, 0.0));

        // Read list.dat
        std::vector<PathListEntry> path_list;
        std::vector<std::string> head;
        int ntotal = 0;
        {
            std::ifstream flist(files.list_file);
            if (flist.is_open()) {
                std::string line;
                int nhead = 0;
                while (std::getline(flist, line) && nhead < nheadx) {
                    if (line.find("---") != std::string::npos) break;
                    head.push_back(line);
                    nhead++;
                }
                if (std::getline(flist, line)) { /* skip label */ }
                PathListEntry entry;
                while (flist >> entry.ip >> entry.sig2u) {
                    path_list.push_back(entry);
                }
                ntotal = static_cast<int>(path_list.size());
            }
        }

        // Read phase.pad for ne3
        int ne3 = 0;
        {
            std::ifstream ppad("phase.pad");
            if (ppad.is_open()) {
                int dummy;
                ppad >> dummy >> dummy >> dummy >> ne3;
            }
        }

        // Read feff.pad and xsect.bin
        FeffPadData pad = read_feff_pad(files.pad_file);
        XsectData xs = read_xsect_bin(p.s02, p.mbconv);

        // Combined title
        int ntitle = xs.ntitle;
        if (ntfms == 1 && ntitle < nheadx)
            xs.title[ntitle++] = titfms;
        for (const auto& h : head)
            if (ntitle < nheadx) xs.title[ntitle++] = h;

        // Write feffnnnn.dat
        if (p.ipr6 == 3) {
            std::vector<int> iplst(ntotal);
            for (int i = 0; i < ntotal; ++i) iplst[i] = path_list[i].ip;
            feffdt(ntotal, iplst.data(), pad.nptot, ntitle, xs.title.data(),
                   pad.ne, pad.iorder, pad.ilinit, pad.rnrmav, pad.edge,
                   pad.potlbl.data(), pad.iz.data(), pad.phc.data(),
                   pad.ck.data(), pad.xk.data(), pad.index.data(),
                   pad.nleg.data(), pad.deg.data(), pad.reff.data(),
                   pad.crit.data(), pad.ipot, pad.rat, pad.achi, pad.phchi);
        }

        // Apply vicorr
        if (std::abs(p.vicorr) >= eps4) {
            for (int ie = 0; ie < pad.ne; ++ie) {
                FeffComplex ck_d(pad.ck[ie].real(), pad.ck[ie].imag());
                FeffComplex ckp = std::sqrt(ck_d * ck_d + 2.0 * coni * p.vicorr);
                double xlam0 = pad.ck[ie].imag() - ckp.imag();
                for (int ipath = 0; ipath < pad.nptot; ++ipath) {
                    pad.achi[ipath][ie] *= static_cast<float>(
                        std::exp(2.0 * pad.reff[ipath] * xlam0));
                }
            }
        }

        // vrcorr edge shift
        float edge = pad.edge;
        if (std::abs(p.vrcorr) > eps4) edge -= static_cast<float>(p.vrcorr);

        // Build k' grid
        double xkp[nex];
        for (int i = 0; i < pad.ne; ++i) {
            double temp = pad.xk[i] * std::abs(static_cast<double>(pad.xk[i])) +
                          2.0 * p.vrcorr;
            xkp[i] = (temp >= 0.0) ? std::sqrt(temp) : -std::sqrt(-temp);
        }

        bool dwcorr = (p.tk > 1.0e-3);
        int ne = pad.ne;
        int ne1 = xs.ne1;

        // Initialize cchi with FMS contribution
        FeffComplex cchi[nex]{};
        for (int ik = 0; ik < ne; ++ik)
            cchi[ik] = p.s02 * gtr[ik];

        // Vicorr broadening
        if (p.vicorr > eps4)
            math::conv(xs.omega.data(), cchi, ne1, p.vicorr);

        // Add DW factors (ispec forced to 3 for anomalous)
        FF2xParams p_afs = p;
        p_afs.ispec = 3;
        int nused = 0;
        dwadd(ntotal, pad.nptot, p_afs, path_list, pad, xs,
              ne, xkp, xkp, cchi, iabs, nused);

        // Configuration average
        FeffComplex chia[nex]{};
        for (int ik = 0; ik < ne; ++ik)
            chia[ik] += cchi[ik] / static_cast<double>(p.nabs);

        if (iabs == p.nabs) {
            std::ofstream xmu_out(files.xmu_file);
            xmu_out << "# " << nused << "/" << ntotal << " paths used\n";

            double rchtot[nex];
            for (int ik = 0; ik < ne; ++ik)
                rchtot[ik] = chia[ik].imag();

            double efermi = edge + xs.omega[0] - xs.emxs[0].real();

            // Convolution
            if (p.mbconv > 0) {
                double wp_half = xs.wp / 2.0;
                exconv(xs.omega.data(), ne1, efermi, xs.s02p, xs.erelax,
                       wp_half, xs.xsnorm.data());
                exconv(xs.omega.data(), ne1, efermi, xs.s02p, xs.erelax,
                       wp_half, rchtot);
            }

            double edg50 = efermi + 50.0 / hart;
            double xsedge = 0.0;
            math::terp(xs.omega.data(), xs.xsnorm.data(), ne1, 1, edg50, xsedge);
            if (p.absolu == 1) xsedge = 1.0;

            xmu_out << "# xsedge+ 50, used to normalize mu "
                    << std::scientific << xsedge << "\n";
            xmu_out << "# " << std::string(71, '-') << "\n";
            xmu_out << "# omega    e    k    mu    mu0     chi     @#\n";

            // Transform cross section to f" units
            for (int ie = 0; ie < ne; ++ie) {
                double energy = xs.emxs[ie].real() + efermi;
                double prefac = 4.0 * pi * alpinv / energy * bohr * bohr;
                xs.xsec[ie] = xs.xsec[ie] / prefac * alpinv * alpinv;
                xs.xsnorm[ie] = xs.xsnorm[ie] / prefac * alpinv * alpinv;
            }

            FeffComplex kxsec[nex];
            for (int i = 0; i < nex; ++i)
                kxsec[i] = cross ? FeffComplex(0.0, 0.0) : xs.xsec[i];

            // f' correction using fprime
            int ne2 = ne - ne1 - ne3;
            FeffComplex cchi_corr[nex]{};
            double fpp[nex];

            fprime(efermi, xs.emxs.data(), ne1, ne3, ne, xs.ik0,
                   kxsec, xs.xsnorm.data(), chia,
                   p.vrcorr, p.vicorr, cchi_corr);

            for (int ie = 0; ie < ne1; ++ie) {
                fpp[ie] = xs.xsnorm[ie] + (xs.xsnorm[ie] * chia[ie]).imag();
                rchtot[ie] = (xs.xsnorm[ie] * chia[ie] + cchi_corr[ie]).real();
            }

            // Second pass for bare f'
            for (int ie = 0; ie < ne; ++ie) chia[ie] = FeffComplex(0.0, 0.0);
            fprime(efermi, xs.emxs.data(), ne1, ne3, ne, xs.ik0,
                   kxsec, xs.xsnorm.data(), chia,
                   p.vrcorr, p.vicorr, cchi_corr);

            // Write output
            for (int ie = 0; ie < ne1; ++ie) {
                double em0 = xs.emxs[ie].real();
                double xsec0 = cchi_corr[ie].real();
                double chi0 = rchtot[ie] - xsec0;

                char buf[128];
                if (ne2 > 0) {
                    // DANES: signs comply with Cromer-Liberman notation
                    std::snprintf(buf, sizeof(buf),
                                  " %10.3f%11.3f%8.3f%13.5e%13.5e%13.5e",
                                  xs.omega[ie] * hart, em0 * hart, xkp[ie] / bohr,
                                  -rchtot[ie], -xsec0, -chi0);
                } else {
                    // FPRIME
                    std::snprintf(buf, sizeof(buf),
                                  " %10.3f%11.3f%13.5e%13.5e%13.5e%13.5e",
                                  xs.omega[ie] * hart, em0 * hart,
                                  -rchtot[ie], -xsec0, fpp[ie], xs.xsnorm[ie]);
                }
                xmu_out << buf << "\n";
            }

            xmu_out.close();
        }
    }
}

} // namespace feff::ff2x
