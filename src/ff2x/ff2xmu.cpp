// XANES mu(E) calculation using FMS+Paths method.
// Converted from: src/FF2X/ff2xmu.f
// Uses FMS Green's function plus path contributions with Debye-Waller factors.

#include "ff2xmu.hpp"
#include "ff2x.hpp"
#include "rdfbin.hpp"
#include "rdxbin.hpp"
#include "feffdt.hpp"
#include "exconv.hpp"
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

void ff2xmu(const FF2xParams& p, int iabs) {
    auto& log = common::logger();

    // Read FMS Green's function from fms.bin
    int nip = p.ipmax - p.ipmin + 1;
    std::vector<std::vector<FeffComplex>> gtrful(nip);
    int ne_fms = 0, ne1_fms = 0, ne3_fms = 0;
    int ntfms = 0;
    std::string titfms;

    {
        std::ifstream fms_in("fms.bin");
        if (fms_in.is_open()) {
            ntfms = 1;
            std::getline(fms_in, titfms);

            int nph_fms = 0, npadx_fms = 0;
            {
                std::string line;
                std::getline(fms_in, line);
                std::istringstream iss(line);
                iss >> ne_fms >> ne1_fms >> ne3_fms >> nph_fms >> npadx_fms;
            }

            // Read all gtr data as one flat complex array
            int total_pts = ne_fms * nip;
            std::vector<FeffComplex> gtrtemp(total_pts);
            common::read_pad_complex(fms_in, npadx_fms,
                                     gtrtemp.data(), total_pts);

            // Distribute into per-spectrum arrays
            int idx = 0;
            for (int iip = 0; iip < nip; ++iip) {
                gtrful[iip].resize(ne_fms);
                for (int j = 0; j < ne_fms; ++j) {
                    gtrful[iip][j] = gtrtemp[idx + j];
                }
                idx += ne_fms;
            }
        }
    }

    // Read xsect.bin
    XsectData xs = read_xsect_bin(p.s02, p.mbconv);

    // Loop over spectrum indices
    for (int iip = p.ipmin; iip <= p.ipmax; iip += p.ipstep) {
        bool cross = is_cross_spectrum(iip);
        auto files = make_spectrum_files(iip);

        // Get gtr for this spectrum
        int iip_idx = iip - p.ipmin;
        std::vector<FeffComplex> gtr;
        if (ntfms > 0 && iip_idx < nip) {
            gtr = gtrful[iip_idx];
        } else {
            gtr.resize(nex, FeffComplex(0.0, 0.0));
        }

        // Read list.dat
        std::vector<PathListEntry> path_list;
        std::vector<std::string> head;
        int ntotal = 0;
        {
            std::ifstream flist(files.list_file);
            if (flist.is_open()) {
                std::string line;
                // Read header
                int nhead = 0;
                while (std::getline(flist, line) && nhead < nheadx) {
                    if (line.find("---") != std::string::npos ||
                        (line.size() > 0 && std::isdigit(line[0]))) break;
                    head.push_back(line);
                    nhead++;
                }
                // Skip label
                if (std::getline(flist, line)) { /* skip */ }
                PathListEntry entry;
                while (flist >> entry.ip >> entry.sig2u) {
                    path_list.push_back(entry);
                }
                ntotal = static_cast<int>(path_list.size());
            }
        }

        // Read feff.pad
        FeffPadData pad = read_feff_pad(files.pad_file);

        // Make combined title
        int ntitle = xs.ntitle;
        if (ntfms == 1 && ntitle < nheadx) {
            xs.title[ntitle++] = titfms;
        }
        for (const auto& h : head) {
            if (ntitle < nheadx)
                xs.title[ntitle++] = h;
        }

        // Write feffnnnn.dat if requested
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
                FeffComplex ckp = std::sqrt(ck_d * ck_d + coni * 2.0 * p.vicorr);
                double xlam0 = pad.ck[ie].imag() - ckp.imag();
                for (int ipath = 0; ipath < pad.nptot; ++ipath) {
                    pad.achi[ipath][ie] *= static_cast<float>(
                        std::exp(2.0 * pad.reff[ipath] * xlam0));
                }
            }
        }

        // Apply vrcorr edge shift
        float edge = pad.edge;
        if (std::abs(p.vrcorr) > eps4) {
            edge -= static_cast<float>(p.vrcorr);
        }

        // Build k' grid (on original energy grid, not fine grid)
        double xkp[nex];
        for (int i = 0; i < pad.ne; ++i) {
            double temp = pad.xk[i] * std::abs(static_cast<double>(pad.xk[i])) +
                          2.0 * p.vrcorr;
            if (temp >= 0.0)
                xkp[i] = std::sqrt(temp);
            else
                xkp[i] = -std::sqrt(-temp);
        }

        bool dwcorr = (p.tk > 1.0e-3);
        int ne = pad.ne;
        int ne1 = xs.ne1;

        // Initialize cchi with FMS contribution
        FeffComplex cchi[nex]{};
        for (int ik = 0; ik < ne; ++ik) {
            cchi[ik] = p.s02 * gtr[ik];
        }

        // Add DW factors (using ne as nkx, xkp as both grids)
        int nused = 0;
        dwadd(ntotal, pad.nptot, p, path_list, pad, xs,
              ne, xkp, xkp, cchi, iabs, nused);

        // Configuration average
        FeffComplex chia[nex]{};
        for (int ik = 0; ik < ne; ++ik) {
            chia[ik] += cchi[ik] / static_cast<double>(p.nabs);
        }

        if (iabs == p.nabs) {
            // Open output
            std::ofstream xmu_out(files.xmu_file);

            xmu_out << "# " << nused << "/" << ntotal << " paths used\n";

            double rchtot[nex];
            for (int ik = 0; ik < ne; ++ik) {
                rchtot[ik] = chia[ik].imag();
            }

            double efermi = edge + xs.omega[0] - xs.emxs[0].real();

            // Convolution with excitation spectrum
            if (p.mbconv > 0) {
                double wp_half = xs.wp / 2.0;
                exconv(xs.omega.data(), ne1, efermi, xs.s02p, xs.erelax,
                       wp_half, xs.xsnorm.data());
                exconv(xs.omega.data(), ne1, efermi, xs.s02p, xs.erelax,
                       wp_half, rchtot);
            }

            // Normalize to xsec at 50 eV above edge
            double edg50 = efermi + 50.0 / hart;
            if (p.ispec == 2) edg50 = efermi;
            double xsedge = 0.0;
            math::terp(xs.omega.data(), xs.xsnorm.data(), ne1, 1, edg50, xsedge);
            if (p.absolu == 1) xsedge = 1.0;

            xmu_out << "# xsedge+ 50, used to normalize mu "
                    << std::scientific << xsedge << "\n";
            xmu_out << "# " << std::string(71, '-') << "\n";
            xmu_out << "# omega    e    k    mu    mu0     chi     @#\n";

            // Cross-section selection
            FeffComplex kxsec[nex];
            for (int i = 0; i < nex; ++i) {
                kxsec[i] = cross ? FeffComplex(0.0, 0.0) : xs.xsec[i];
            }

            // Brouder method correction
            FeffComplex cchi_corr[nex]{};
            double vi0 = 0.0;
            xscorr(p.ispec, xs.emxs.data(), ne1, ne, xs.ik0,
                   kxsec, xs.xsnorm.data(), chia,
                   p.vrcorr, vi0, cchi_corr);

            for (int ie = 0; ie < ne1; ++ie) {
                rchtot[ie] = (kxsec[ie] + xs.xsnorm[ie] * chia[ie] + cchi_corr[ie]).imag();
            }

            // Second pass with zeroed chia for bare xsec
            for (int ie = 0; ie < ne; ++ie)
                chia[ie] = FeffComplex(0.0, 0.0);
            xscorr(p.ispec, xs.emxs.data(), ne1, ne, xs.ik0,
                   kxsec, xs.xsnorm.data(), chia,
                   p.vrcorr, vi0, cchi_corr);

            for (int ie = 0; ie < ne1; ++ie) {
                cchi_corr[ie] = FeffComplex(rchtot[ie],
                                            (kxsec[ie] + cchi_corr[ie]).imag());
            }

            if (p.vicorr > eps4 && ntotal == 0) {
                math::conv(xs.omega.data(), cchi_corr, ne1, p.vicorr);
            }

            // Write output
            for (int ie = 0; ie < ne1; ++ie) {
                double em0 = xs.emxs[ie].real();
                double xsec0 = cchi_corr[ie].imag();
                rchtot[ie] = cchi_corr[ie].real();
                double chi0 = (rchtot[ie] - xsec0) / xsedge;

                char buf[128];
                std::snprintf(buf, sizeof(buf), " %10.3f%11.3f%8.3f%13.5e%13.5e%13.5e",
                              xs.omega[ie] * hart, em0 * hart, xkp[ie] / bohr,
                              rchtot[ie] / xsedge, xsec0 / xsedge, chi0);
                xmu_out << buf << "\n";
            }

            xmu_out.close();
        }
    }
}

} // namespace feff::ff2x
