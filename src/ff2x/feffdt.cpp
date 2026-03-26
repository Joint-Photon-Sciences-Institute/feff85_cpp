// Write feffNNNN.dat and files.dat output files.
// Converted from: src/FF2X/feffdt.f, fdthea.f, fdtarr.f

#include "feffdt.hpp"

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <feff/feff_input.hpp>

#include "../common/logging.hpp"

#include "../math/phase_amplitude.hpp"

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace feff::ff2x {

static constexpr double eps_tiny = 1.0e-16;

// Pi-jump correction: remove 2*pi jumps in phase
static void pijump(double& phase, double phase_old) {
    while (phase - phase_old > pi) phase -= 2.0 * pi;
    while (phase - phase_old < -pi) phase += 2.0 * pi;
}

void fdtarr(int ne, float reff, int lzero,
            const float achi[], const float phchi[],
            const std::complex<float> caps[], const float xk[],
            const std::complex<float> ck[],
            double col1[], double col2[], double col3[], double col4[],
            double col5[], double col6[], double col7[]) {

    double phffo = 0.0, cdelto = 0.0;

    for (int ie = 0; ie < ne; ++ie) {
        // cchi = achi * exp(i * phchi)
        FeffComplex cchi_val = static_cast<double>(achi[ie]) *
                               std::exp(coni * static_cast<double>(phchi[ie]));

        double xlam = 1.0e10;
        double ck_imag = static_cast<double>(ck[ie].imag());
        if (std::abs(ck_imag) > eps_tiny)
            xlam = 1.0 / ck_imag;

        double redfac = std::exp(-2.0 * static_cast<double>(caps[ie].imag()));
        double cdelt = 2.0 * static_cast<double>(caps[ie].real());

        double reff_d = static_cast<double>(reff);
        double xk_d = static_cast<double>(xk[ie]);

        FeffComplex cfms = cchi_val * xk_d * reff_d * reff_d *
                           std::exp(2.0 * reff_d / xlam) / redfac;

        double phff = 0.0;
        if (std::abs(cchi_val) >= eps_tiny) {
            phff = std::atan2(cchi_val.imag(), cchi_val.real());
        }

        // Remove 2*pi jumps in phases
        if (ie > 0) {
            pijump(phff, phffo);
            pijump(cdelt, cdelto);
        }
        phffo = phff;
        cdelto = cdelt;

        col1[ie] = xk_d / bohr;
        col2[ie] = cdelt + lzero * pi;
        col3[ie] = std::abs(cfms) * bohr;
        col4[ie] = phff - cdelt - lzero * pi;
        col5[ie] = redfac;
        col6[ie] = xlam * bohr;
        col7[ie] = static_cast<double>(ck[ie].real()) / bohr;
    }
}

void fdthea(int ntext, const std::string text[], int ip, int iorder,
            int nleg_val, float deg, float reff, float rnrmav, float edge,
            const std::array<float, 3> rat_leg[], const int ipot_leg[],
            const int iz[], const std::string potlbl[],
            int& nlines, std::vector<std::string>& lines) {

    lines.clear();
    char buf[256];

    // Title lines
    for (int i = 0; i < ntext; ++i) {
        lines.push_back(" " + text[i]);
    }

    // Path line
    std::snprintf(buf, sizeof(buf), " Path%5d      icalc %7d", ip, iorder);
    lines.push_back(buf);

    // Separator
    lines.push_back(" " + std::string(71, '-'));

    // nleg, deg, reff, rnrmav, edge
    std::snprintf(buf, sizeof(buf), " %3d%8.3f%9.4f%10.4f%11.5f nleg, deg, reff, rnrmav(bohr), edge",
                  nleg_val, static_cast<double>(deg),
                  static_cast<double>(reff) * bohr,
                  static_cast<double>(rnrmav),
                  static_cast<double>(edge) * hart);
    lines.push_back(buf);

    // Column headers
    lines.push_back("        x         y         z   pot at#");

    // Absorbing atom (last leg = index nleg-1 in 0-indexed)
    int iabs = nleg_val - 1;
    std::snprintf(buf, sizeof(buf), " %9.4f%10.4f%10.4f%3d%4d %s   absorbing atom",
                  rat_leg[iabs][0] * bohr,
                  rat_leg[iabs][1] * bohr,
                  rat_leg[iabs][2] * bohr,
                  ipot_leg[iabs],
                  iz[ipot_leg[iabs]],
                  potlbl[ipot_leg[iabs]].c_str());
    lines.push_back(buf);

    // Scattering atoms
    for (int ileg = 0; ileg < nleg_val - 1; ++ileg) {
        std::snprintf(buf, sizeof(buf), " %9.4f%10.4f%10.4f%3d%4d %s",
                      rat_leg[ileg][0] * bohr,
                      rat_leg[ileg][1] * bohr,
                      rat_leg[ileg][2] * bohr,
                      ipot_leg[ileg],
                      iz[ipot_leg[ileg]],
                      potlbl[ipot_leg[ileg]].c_str());
        lines.push_back(buf);
    }

    // Column labels
    lines.push_back("    k   real[2*phc]   mag[feff]  phase[feff]"
                     " red factor   lambda     real[p]@#");

    nlines = static_cast<int>(lines.size());
}

void feffdt(int ntotal, const int iplst[], int nptot,
            int ntext, const std::string text[],
            int ne, int iorder, int ilinit, float rnrmav, float edge,
            const std::string potlbl[], const int iz[],
            const std::complex<float> phc[], const std::complex<float> ck[],
            const float xk[], const int index[], const int nleg[],
            const float deg[], const float reff[], const float crit[],
            const std::vector<std::vector<int>>& ipot,
            const std::vector<std::vector<std::array<float, 3>>>& rat,
            const std::vector<std::vector<float>>& achi,
            const std::vector<std::vector<float>>& phchi) {

    auto& log = common::logger();

    log.wlog(" feffdt, feff.pad to feff.dat conversion: " +
             vfeff_default + "release " + vf85e_default);

    for (int i = 0; i < ntext; ++i) {
        log.wlog(" " + text[i]);
    }

    {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%8d paths to process", nptot);
        log.wlog(buf);
    }

    // Open files.dat
    std::ofstream files_out("files.dat");
    if (!files_out.is_open()) {
        throw std::runtime_error("Cannot open files.dat for writing");
    }

    for (int i = 0; i < ntext; ++i) {
        files_out << " " << text[i] << "\n";
    }
    files_out << " " << std::string(71, '-') << "\n";
    files_out << "    file        sig2   amp ratio    "
              << "deg    nlegs  r effective\n";

    log.wlog("    path     filename");

    for (int ilist = 0; ilist < ntotal; ++ilist) {
        // Find path index
        int i = -1;
        for (int j = 0; j < nptot; ++j) {
            if (iplst[ilist] == index[j]) {
                i = j;
                break;
            }
        }
        if (i < 0) {
            char buf[128];
            std::snprintf(buf, sizeof(buf), " did not find path ilist=%d, iplst=%d",
                          ilist, iplst[ilist]);
            log.wlog(buf);
            continue;
        }

        // Prepare filename
        char fname[16];
        std::snprintf(fname, sizeof(fname), "feff%04d.dat", index[i]);

        {
            char buf[64];
            std::snprintf(buf, sizeof(buf), "%8d     %s", i, fname);
            log.wlog(buf);
        }

        // Write to files.dat
        {
            char buf[128];
            std::snprintf(buf, sizeof(buf), " %s%8.5f%10.3f%10.3f%6d%9.4f",
                          fname, 0.0, static_cast<double>(crit[i]),
                          static_cast<double>(deg[i]),
                          nleg[i], static_cast<double>(reff[i]) * bohr);
            files_out << buf << "\n";
        }

        // Build header
        int nlines = 0;
        std::vector<std::string> lines;
        fdthea(ntext, text, index[i], iorder, nleg[i], deg[i], reff[i],
               rnrmav, edge, rat[i].data(), ipot[i].data(), iz, potlbl,
               nlines, lines);

        // Compute columns
        double col1[nex], col2[nex], col3[nex], col4[nex], col5[nex], col6[nex], col7[nex];
        fdtarr(ne, reff[i], ilinit, achi[i].data(), phchi[i].data(),
               phc, xk, ck, col1, col2, col3, col4, col5, col6, col7);

        // Write feffNNNN.dat
        std::ofstream fout(fname);
        if (!fout.is_open()) {
            throw std::runtime_error(std::string("Cannot open ") + fname);
        }

        for (const auto& line : lines) {
            fout << line << "\n";
        }

        fout << std::fixed;
        for (int ie = 0; ie < ne; ++ie) {
            char dbuf[128];
            std::snprintf(dbuf, sizeof(dbuf),
                          " %6.3f %11.4e %11.4e %11.4e %10.3e %11.4e %11.4e",
                          col1[ie], col2[ie], col3[ie], col4[ie],
                          col5[ie], col6[ie], col7[ie]);
            fout << dbuf << "\n";
        }
        fout.close();
    }

    files_out.close();
}

} // namespace feff::ff2x
