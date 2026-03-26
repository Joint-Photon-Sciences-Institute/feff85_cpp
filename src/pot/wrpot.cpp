// Write pot.pad output file in PAD (Packed ASCII Data) format.
// Converted from: src/POT/wrpot.f
//
// Opens pot.pad and writes all potential/density data needed by
// downstream modules (XSPH, etc.).

#include "wrpot.hpp"

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include "../common/pad_io.hpp"
#include "../common/logging.hpp"
#include "../common/string_utils.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstring>

namespace feff::pot {

void wrpot(int nph, int ntitle, const char title[][80],
           double rnrmav, double xmu, double vint, double rhoint,
           double emu, double s02, double erelax, double wp,
           double ecv, double rs, double xf, double qtotel,
           const int* imt, const double* rmt, const int* inrm, const double* rnrm,
           const double* folp, const double* folpx, const double* xnatph,
           const double* dgc0, const double* dpc0,
           const double* dgc, const double* dpc,
           const double* adgc, const double* adpc,
           const double* edens, const double* vclap, const double* vtot,
           const double* edenvl, const double* vvalgs, const double* dmag,
           const double* xnval,
           const double* eorb, const int* kappa, const int* iorb,
           const double* qnrm, const double* xnmues,
           int nohole, int ihole,
           int inters, double totvol, int iafolp,
           const double* xion, int iunf, const int* iz, int jumprm)
{
    constexpr int npadx = 8;

    std::ofstream out("pot.pad", std::ios::out);
    feff::common::check_file_open(out, "pot.pad", "wrpot");

    // Write header line: ntitle, nph, npadx, nohole, ihole, inters, iafolp, jumprm, iunf
    out << std::setw(5) << ntitle
        << std::setw(5) << nph
        << std::setw(5) << npadx
        << std::setw(5) << nohole
        << std::setw(5) << ihole
        << std::setw(5) << inters
        << std::setw(5) << iafolp
        << std::setw(5) << jumprm
        << std::setw(5) << iunf
        << "\n";

    // Write title lines
    for (int i = 0; i < ntitle; i++) {
        // Find effective length (trim trailing spaces)
        int ll = 79;
        while (ll >= 0 && (title[i][ll] == ' ' || title[i][ll] == '\0')) ll--;
        ll++;
        out.write(title[i], ll);
        out << "\n";
    }

    // Misc scalar data packed as doubles
    double dum[13];
    dum[0]  = rnrmav;
    dum[1]  = xmu;
    dum[2]  = vint;
    dum[3]  = rhoint;
    dum[4]  = emu;
    dum[5]  = s02;
    dum[6]  = erelax;
    dum[7]  = wp;
    dum[8]  = ecv;
    dum[9]  = rs;
    dum[10] = xf;
    dum[11] = qtotel;
    dum[12] = totvol;
    feff::common::write_pad_double(out, npadx, dum, 13);

    // imt array (0:nph)
    for (int i = 0; i <= nph; i++) {
        out << std::setw(5) << imt[i];
    }
    out << "\n";

    // rmt array
    feff::common::write_pad_double(out, npadx, rmt, nph + 1);

    // inrm array
    for (int i = 0; i <= nph; i++) {
        out << std::setw(5) << inrm[i];
    }
    out << "\n";

    // iz array
    for (int i = 0; i <= nph; i++) {
        out << std::setw(5) << iz[i];
    }
    out << "\n";

    // kappa array (1:30) - first 30 elements from absorber
    for (int i = 0; i < 30; i++) {
        out << std::setw(5) << kappa[i];
    }
    out << "\n";

    // Double arrays via PAD
    feff::common::write_pad_double(out, npadx, rnrm, nph + 1);
    feff::common::write_pad_double(out, npadx, folp, nph + 1);
    feff::common::write_pad_double(out, npadx, folpx, nph + 1);
    feff::common::write_pad_double(out, npadx, xnatph, nph + 1);
    feff::common::write_pad_double(out, npadx, xion, nph + 1);
    feff::common::write_pad_double(out, npadx, dgc0, 251);
    feff::common::write_pad_double(out, npadx, dpc0, 251);
    feff::common::write_pad_double(out, npadx, dgc, 251 * 30 * (nph + 1));
    feff::common::write_pad_double(out, npadx, dpc, 251 * 30 * (nph + 1));
    feff::common::write_pad_double(out, npadx, adgc, 10 * 30 * (nph + 1));
    feff::common::write_pad_double(out, npadx, adpc, 10 * 30 * (nph + 1));
    feff::common::write_pad_double(out, npadx, edens, 251 * (nph + 1));
    feff::common::write_pad_double(out, npadx, vclap, 251 * (nph + 1));
    feff::common::write_pad_double(out, npadx, vtot, 251 * (nph + 1));
    feff::common::write_pad_double(out, npadx, edenvl, 251 * (nph + 1));
    feff::common::write_pad_double(out, npadx, vvalgs, 251 * (nph + 1));
    feff::common::write_pad_double(out, npadx, dmag, 251 * (nph + 1));
    feff::common::write_pad_double(out, npadx, xnval, 30 * (nph + 1));
    feff::common::write_pad_double(out, npadx, eorb, 30);

    // iorb array: (-4:3, 0:nphx) -> 8 ints per iph
    for (int iph = 0; iph <= nph; iph++) {
        for (int i = 0; i < 8; i++) {
            // iorb is stored as iorb[8*iph + i] where i=0 corresponds to kappa=-4
            out << std::setw(3) << iorb[i + 8 * iph];
        }
        out << "\n";
    }

    feff::common::write_pad_double(out, npadx, qnrm, nph + 1);
    feff::common::write_pad_double(out, npadx, xnmues, (lx + 1) * (nph + 1));

    out.close();
}

} // namespace feff::pot
