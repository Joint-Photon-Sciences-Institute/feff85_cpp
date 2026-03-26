// Automatic overlap fractions.
// Converted from src/POT/afolp.f

#include "afolp.hpp"
#include "istprm.hpp"
#include "../common/logging.hpp"
#include <feff/constants.hpp>
#include <sstream>
#include <iomanip>

namespace feff::pot {

void afolp(bool verbse, int nph, int nat, const int* iphat, const double* rat,
           const int* iatph, const double* xnatph,
           const int* novr, const int* iphovr, const int* nnovr, const double* rovr,
           double* folp, double* folpx, int iafolp,
           double* edens, double* edenvl,
           double* dmag, double* vclap, double* vtot, double* vvalgs,
           int* imt, int* inrm, double* rmt, double* rnrm,
           int ixc, double& rhoint, double& vint, double& rs, double& xf,
           double xmu, double& xmunew,
           double& rnrmav, double& qtotel, int& inters, double totvol)
{
    // Save un-folped rmt values
    double rmtx[nphx + 1];
    for (int iph = 0; iph <= nph; ++iph) {
        rmtx[iph] = rmt[iph] / folp[iph];
    }

    if (verbse) {
        feff::common::logger().wlog(
            "    : ipot, Norman radius, Muffin tin radius, Overlap");
    }

    if (iafolp >= 0) {
        for (int iph = 0; iph <= nph; ++iph) {
            folp[iph] = folpx[iph];
            rmt[iph] = folp[iph] * rmtx[iph];

            if (verbse) {
                std::ostringstream slog;
                slog << std::setw(10) << iph
                     << std::scientific << std::setprecision(5)
                     << std::setw(13) << rnrm[iph] * bohr
                     << std::setw(13) << rmt[iph] * bohr
                     << std::setw(13) << folp[iph];
                feff::common::logger().wlog(slog.str());
            }
        }

        int idmag = 0;
        istprm(nph, nat, iphat, rat, iatph, xnatph,
               novr, iphovr, nnovr, rovr, folp, folpx, iafolp,
               edens, edenvl, idmag,
               dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,
               ixc, rhoint, vint, rs, xf, xmu, xmunew,
               rnrmav, qtotel, inters, totvol);
    }
}

} // namespace feff::pot
