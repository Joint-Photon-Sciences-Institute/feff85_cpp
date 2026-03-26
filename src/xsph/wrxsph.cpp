// Write phase.pad output in PAD format.
// Converted from src/XSPH/wrxsph.f

#include "wrxsph.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include "../common/pad_io.hpp"
#include <fstream>
#include <iomanip>
#include <vector>

namespace feff::xsph {

void wrxsph(const std::string& phpad,
            int nsp, int ne, int ne1, int ne3, int nph, int ihole,
            double rnrmav, double xmu, double edge,
            int ik0, int ixc, double rs, double vint,
            const FeffComplex em[], const FeffComplex* eref,
            const int lmax[], const int iz[], const char potlbl[][7],
            const FeffComplex* ph, const FeffComplex* rkk) {

    constexpr int npadx = 8;

    std::ofstream out(phpad);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open " + phpad);
    }

    // Write header line
    out << std::setw(5) << nsp
        << std::setw(5) << ne
        << std::setw(5) << ne1
        << std::setw(5) << ne3
        << std::setw(5) << nph
        << std::setw(5) << ihole
        << std::setw(5) << ik0
        << std::setw(5) << npadx
        << std::setw(5) << ixc
        << std::fixed << std::setprecision(5)
        << std::setw(11) << rs
        << std::setw(11) << vint
        << "\n";

    // Write rnrmav, xmu, edge
    double dum[3] = {rnrmav, xmu, edge};
    feff::common::write_pad_double(out, npadx, dum, 3);

    // Write em
    feff::common::write_pad_complex(out, npadx, em, ne);

    // Write eref
    int ii = 0;
    std::vector<FeffComplex> temp(nex * (2 * ltot + 1));
    for (auto& t : temp) t = FeffComplex(0.0, 0.0);

    ii = 0;
    for (int isp = 0; isp < nsp; isp++) {
        for (int ie = 0; ie < ne; ie++) {
            temp[ii++] = eref[ie + nex * isp];
        }
    }
    feff::common::write_pad_complex(out, npadx, temp.data(), ii);

    // Write phase shifts for each potential
    for (int iph = 0; iph <= nph; iph++) {
        out << std::setw(4) << lmax[iph]
            << std::setw(4) << iz[iph]
            << " " << std::string(potlbl[iph], 6) << "\n";

        for (int isp = 0; isp < nsp; isp++) {
            ii = 0;
            for (int ie = 0; ie < ne; ie++) {
                for (int ll = -lmax[iph]; ll <= lmax[iph]; ll++) {
                    temp[ii++] = ph[ie + nex * ((ll + ltot) + (2 * ltot + 1) * (isp + nspx * iph))];
                }
            }
            feff::common::write_pad_complex(out, npadx, temp.data(), ii);
        }
    }

    // Write rkk
    // rkk is stored in C row-major order as rkk[isp * nex * 8 + ie * 8 + kdif]
    // (filled by xsect which receives FeffComplex rkk[][8]).
    // We write in Fortran column-major order: inner ie, then kdif, then isp.
    ii = 0;
    for (int isp = 0; isp < nsp; isp++) {
        for (int kdif = 0; kdif < 8; kdif++) {
            for (int ie = 0; ie < ne; ie++) {
                temp[ii++] = rkk[isp * nex * 8 + ie * 8 + kdif];
            }
        }
    }
    feff::common::write_pad_complex(out, npadx, temp.data(), ii);
}

} // namespace feff::xsph
