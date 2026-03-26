// EXAFS chi(k) calculation using MS paths expansion.
// Converted from: src/FF2X/ff2chi.f
// Adds contributions from each path and absorber, including
// Debye-Waller factors. Writes chi.dat and xmu.dat.

#include "ff2chi.hpp"
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

#include "../math/interpolation.hpp"

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace feff::ff2x {

// getxk: compute signed k from energy
static double getxk(double e) {
    if (e >= 0.0) return std::sqrt(2.0 * e);
    return -std::sqrt(-2.0 * e);
}

void ff2chi(const FF2xParams& p, int iabs) {
    auto& log = common::logger();

    // Check for phase.pad to determine ne3 for DANES skip
    int ne3 = 0;
    {
        std::ifstream ppad("phase.pad");
        if (ppad.is_open() && std::abs(p.ispec) == 3) {
            int dummy;
            ppad >> dummy >> dummy >> dummy >> ne3;
        }
    }

    // Loop over spectrum indices (ipmin to ipmax)
    for (int iip = p.ipmin; iip <= p.ipmax; iip += p.ipstep) {
        bool cross = is_cross_spectrum(iip);
        auto files = make_spectrum_files(iip);

        // Read list.dat
        std::vector<PathListEntry> path_list;
        std::vector<std::string> head;
        {
            std::ifstream flist(files.list_file);
            if (!flist.is_open()) {
                throw std::runtime_error("Cannot open " + files.list_file);
            }
            // Read header
            std::string line;
            while (std::getline(flist, line)) {
                if (line.empty() || line[0] != ' ') break;
                head.push_back(line);
                if (line.find("---") != std::string::npos) break;
            }
            // Skip label line
            if (std::getline(flist, line)) { /* skip */ }
            // Read path list (one line per path, first two fields only)
            while (std::getline(flist, line)) {
                if (line.empty()) continue;
                std::istringstream iss(line);
                PathListEntry entry;
                if (iss >> entry.ip >> entry.sig2u) {
                    path_list.push_back(entry);
                }
            }
        }
        int ntotal = static_cast<int>(path_list.size());

        // Read xsect.bin
        XsectData xs = read_xsect_bin(p.s02, p.mbconv);

        // Read feff.pad
        FeffPadData pad = read_feff_pad(files.pad_file);

        // Make combined title
        int ntitle = xs.ntitle;
        for (const auto& h : head) {
            if (ntitle < nheadx) {
                xs.title[ntitle++] = h;
            }
        }

        // Write feffnnnn.dat if requested
        if (p.ipr6 == 3) {
            std::vector<int> iplst(ntotal);
            for (int i = 0; i < ntotal; ++i) iplst[i] = path_list[i].ip;
            feffdt(ntotal, iplst.data(), pad.nptot, ntitle, xs.title.data(),
                   (std::abs(p.ispec) == 3) ? pad.ne : xs.ne1,
                   pad.iorder, pad.ilinit, pad.rnrmav, pad.edge,
                   pad.potlbl.data(), pad.iz.data(), pad.phc.data(),
                   pad.ck.data(), pad.xk.data(), pad.index.data(),
                   pad.nleg.data(), pad.deg.data(), pad.reff.data(),
                   pad.crit.data(), pad.ipot, pad.rat, pad.achi, pad.phchi);
        }

        // Compare grids (xsect vs feff.pad)
        if (iabs == 1) {
            for (int i = 0; i < xs.nxsec; ++i) {
                double del = static_cast<double>(pad.xk[i]) * pad.xk[i] -
                             xs.xkxs[i] * xs.xkxs[i];
                if (std::abs(p.ispec) != 3 && std::abs(del) > 10 * eps4) {
                    log.wlog(" Emesh in feff.pad and xsect.bin different.");
                    log.wlog(" Results may be meaningless, check input files.");
                    throw std::runtime_error("FF2CHI-1: grid mismatch");
                }
            }
        }

        // Apply vicorr mean free path correction
        if (std::abs(p.vicorr) >= eps4) {
            for (int ipath = 0; ipath < pad.nptot; ++ipath) {
                for (int ie = 0; ie < pad.ne; ++ie) {
                    FeffComplex ck_d(pad.ck[ie].real(), pad.ck[ie].imag());
                    FeffComplex ckp = std::sqrt(ck_d * ck_d + coni * 2.0 * p.vicorr);
                    double xlam0 = pad.ck[ie].imag() - ckp.imag();
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

        // Build fine k-grids
        double delk = 0.05 * bohr;
        double tmp = (pad.xk[0] >= 0) ? 1.0 : -1.0;
        double e = tmp * pad.xk[0] * pad.xk[0] / 2.0 + p.vrcorr;
        double xkpmin = getxk(e);
        int n = static_cast<int>(xkpmin / delk);
        if (xkpmin > 0) n++;
        double xkmin = n * delk;

        int ik0p = 0;
        int ne1 = xs.ne1;
        int nkx = 0;

        double xkp[nfinex], xk0[nfinex];
        int ik0 = (std::abs(p.ispec) != 3) ? 0 : xs.ik0;

        for (int i = 0; i < nfinex; ++i) {
            xkp[i] = xkmin + delk * i;
            tmp = (xkp[i] >= 0) ? 1.0 : -1.0;
            e = tmp * xkp[i] * xkp[i] / 2.0 - p.vrcorr;
            xk0[i] = getxk(e);
            if (xk0[i] < eps4) ik0p = i;
            if (xk0[i] > pad.xk[ne1 - 1] + eps4) break;
            nkx = i + 1;
        }

        // Initialize chi accumulator
        FeffComplex cchi[nfinex]{};

        // Add Debye-Waller factors and sum paths
        int nused = 0;
        dwadd(ntotal, pad.nptot, p, path_list, pad, xs,
              nkx, xk0, xkp, cchi, iabs, nused);

        // Configuration average
        FeffComplex chia[nfinex]{};
        if (iabs == 1) {
            for (int ie = 0; ie < nfinex; ++ie)
                chia[ie] = FeffComplex(0.0, 0.0);
        }
        // Note: multi-absorber averaging with chia.bin would go here

        for (int ik = 0; ik < nkx; ++ik) {
            chia[ik] += cchi[ik] / static_cast<double>(p.nabs);
        }

        if (iabs == p.nabs) {
            // Open output files
            std::ofstream chi_out(files.chi_file);
            std::ofstream xmu_out(files.xmu_file);

            // Write headers
            // (simplified -- full header writing omitted for brevity)
            chi_out << "# " << nused << "/" << ntotal << " paths used\n";
            chi_out << "# " << std::string(71, '-') << "\n";
            chi_out << "#      k          chi          mag           phase @#\n";

            // Compute rchtot
            double rchtot[nfinex];
            for (int ik = 0; ik < nkx; ++ik) {
                if (std::abs(p.ispec) != 3)
                    rchtot[ik] = chia[ik].imag();
                else
                    rchtot[ik] = chia[ik].real();
            }

            // Prepare output energy grid
            double efermi = edge + xs.omega[0] - xs.emxs[0].real();
            double omegax[nfinex];
            for (int ik = 0; ik < nkx; ++ik) {
                if (xkp[ik] < 0.0)
                    omegax[ik] = -xkp[ik] * xkp[ik] / 2.0 + efermi;
                else
                    omegax[ik] = xkp[ik] * xkp[ik] / 2.0 + efermi;
            }

            // Convolution with excitation spectrum
            if (p.mbconv > 0) {
                double wp_half = xs.wp / 2.0;
                exconv(xs.omega.data(), ne1, efermi, xs.s02p, xs.erelax,
                       wp_half, xs.xsnorm.data());
                exconv(omegax, nkx, efermi, xs.s02p, xs.erelax,
                       wp_half, rchtot);
            }

            // Write chi.dat
            for (int ik = 0; ik < nkx; ++ik) {
                FeffComplex ccc = chia[ik];
                double phase = 0.0;
                if (std::abs(ccc) > 0.0) {
                    phase = std::atan2(ccc.imag(), ccc.real());
                }
                static double phase0 = 0.0;
                if (ik > 0) {
                    while (phase - phase0 > pi) phase -= 2.0 * pi;
                    while (phase - phase0 < -pi) phase += 2.0 * pi;
                }
                phase0 = phase;

                char buf[128];
                std::snprintf(buf, sizeof(buf), " %10.4f   %13.6e %13.6e %13.6e",
                              xkp[ik] / bohr, rchtot[ik], std::abs(ccc), phase0);
                chi_out << buf << "\n";
            }
            chi_out.close();

            // Write xmu.dat
            double edg50 = efermi + 50.0 / hart;
            double xsedge = 0.0;
            math::terp(xs.omega.data(), xs.xsnorm.data(), ne1, 1, edg50, xsedge);
            if (p.absolu == 1) xsedge = 1.0;

            xmu_out << "# xsedge+50, used to normalize mu " << std::scientific
                    << xsedge << "\n";
            xmu_out << "# " << std::string(71, '-') << "\n";
            xmu_out << "# omega    e    k    mu    mu0     chi     @#\n";

            // Cross-section array selection
            FeffComplex kxsec[nex];
            for (int i = 0; i < nex; ++i) {
                kxsec[i] = cross ? FeffComplex(0.0, 0.0) : xs.xsec[i];
            }

            // Edge correction
            FeffComplex chia_corr[nex]{};
            FeffComplex cchi_corr[nex]{};

            if (std::abs(p.ispec) == 3) {
                // Transform cross section units for anomalous scattering
                for (int ie = 0; ie < pad.ne; ++ie) {
                    double energy = xs.emxs[ie].real() + efermi;
                    double prefac = 4.0 * pi * alpinv / energy * bohr * bohr;
                    kxsec[ie] = kxsec[ie] / prefac * alpinv * alpinv;
                    xs.xsnorm[ie] = xs.xsnorm[ie] / prefac * alpinv * alpinv;
                }
                int ne2 = pad.ne - ne1 - ne3;
                fprime(efermi, xs.emxs.data(), ne1, ne3, pad.ne, xs.ik0,
                       kxsec, xs.xsnorm.data(), chia_corr,
                       p.vrcorr, p.vicorr, cchi_corr);
                for (int ie = 0; ie < ne1; ++ie) {
                    xs.omega[ie] = cchi_corr[ie].real();
                }
            } else {
                xscorr(p.ispec, xs.emxs.data(), ne1, pad.ne, xs.ik0,
                       kxsec, xs.xsnorm.data(), chia_corr,
                       p.vrcorr, p.vicorr, cchi_corr);
                for (int ie = 0; ie < ne1; ++ie) {
                    xs.omega[ie] = (kxsec[ie] + cchi_corr[ie]).imag();
                }
            }

            for (int ik = 0; ik < nkx; ++ik) {
                double em0 = omegax[ik] - efermi + edge;
                double xsec0 = 0.0, xsnor0 = 0.0;
                math::terp(xs.xkxs.data(), xs.omega.data(), ne1, 1, xk0[ik], xsec0);
                math::terp(xs.xkxs.data(), xs.xsnorm.data(), ne1, 1, xk0[ik], xsnor0);

                double chi0;
                if (omegax[ik] >= efermi)
                    chi0 = xsnor0 * rchtot[ik];
                else
                    chi0 = xsnor0 * rchtot[ik0p];

                char buf[128];
                if (std::abs(p.ispec) != 3) {
                    std::snprintf(buf, sizeof(buf), " %10.3f%11.3f%8.3f%13.5e%13.5e%13.5e",
                                  omegax[ik] * hart, em0 * hart, xkp[ik] / bohr,
                                  (chi0 + xsec0) / xsedge,
                                  xsec0 / xsedge, rchtot[ik]);
                } else {
                    std::snprintf(buf, sizeof(buf), " %10.3f%11.3f%8.3f%13.5e%13.5e%13.5e",
                                  omegax[ik] * hart, em0 * hart, xkp[ik] / bohr,
                                  -(xsec0 + chi0), -xsec0, -chi0);
                }
                xmu_out << buf << "\n";
            }
            xmu_out.close();
        }
    }
}

} // namespace feff::ff2x
