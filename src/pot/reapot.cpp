// Read input files for POT module.
// Converted from src/POT/reapot.f

#include "reapot.hpp"
#include <feff/json_io.hpp>
#include <feff/constants.hpp>
#include "../common/string_utils.hpp"
#include <sstream>
#include <iomanip>

namespace feff::pot {

void reapot(int& mpot, double& rgrd, int& ntitle, std::string* title,
            int& ipr1, int& ispec, int& nohole, int& ihole, double& gamach,
            int& nph, int* iz, int* lmaxsc, double* xnatph,
            double* xion, int& iunf, int& ixc, int& jumprm, int& iafolp,
            double* folp, int& inters, double& totvol,
            float& rfms1, int& lfms1, int& nscmt, double& ca1, int& nmix,
            double& ecv, int& icoul,
            int* novr, int* iphovr, int* nnovr, double* rovr,
            int& nat, double* rat, int* iphat, int* iatph)
{
    // Read geometry from geom.json
    int ibounc[natx];
    feff::json_io::read_geom_json(nat, nph, iatph,
                                   reinterpret_cast<double(*)[3]>(rat),
                                   iphat, ibounc);

    // Read potential parameters from pot.json
    feff::json_io::read_pot_json(mpot, nph, ntitle, ihole, ipr1, iafolp,
                                  ixc, ispec, nmix, nohole, jumprm, inters,
                                  nscmt, icoul, lfms1, iunf,
                                  gamach, rgrd, ca1, ecv, totvol, rfms1,
                                  title, iz, lmaxsc, xnatph, xion, folp,
                                  novr,
                                  reinterpret_cast<int(*)[novrx]>(iphovr),
                                  reinterpret_cast<int(*)[novrx]>(nnovr),
                                  reinterpret_cast<double(*)[novrx]>(rovr));

    // Transform to code units (Bohr and Hartrees = atomic units)
    rfms1 = rfms1 / static_cast<float>(bohr);
    gamach = gamach / hart;
    ecv = ecv / hart;
    totvol = totvol / (bohr * bohr * bohr);

    for (int iat = 0; iat < nat; ++iat) {
        for (int i = 0; i < 3; ++i) {
            rat[iat * 3 + i] = rat[iat * 3 + i] / bohr;
        }
    }

    for (int iph = 0; iph <= nph; ++iph) {
        for (int iovr = 0; iovr < novr[iph]; ++iovr) {
            rovr[iph * novrx + iovr] = rovr[iph * novrx + iovr] / bohr;
        }
    }

    // Add lines to the title
    if (mpot == 1) {
        std::string s1, s2, s3;

        if (nat > 1) {
            if (rfms1 < 0) rfms1 = 0;
            if (nscmt > 0) {
                std::ostringstream oss;
                oss << " POT  SCF" << std::setw(4) << nscmt
                    << std::fixed << std::setprecision(4)
                    << std::setw(8) << rfms1 * static_cast<float>(bohr)
                    << std::setw(4) << lfms1;
                s1 = oss.str();
            } else {
                s1 = " POT  Non-SCF";
            }
        } else {
            s1 = " POT  used OVERLAP geometry,";
        }

        if (nohole == 0) {
            s2 = ", NO core-hole,";
        } else if (nohole == 2) {
            s2 = ", screened core-hole,";
        } else {
            s2 = ", core-hole,";
        }

        if (iafolp < 0) {
            std::ostringstream oss;
            oss << " FOLP (folp(0)=" << std::fixed << std::setprecision(3)
                << std::setw(6) << folp[0] << ")";
            s3 = oss.str();
        } else {
            std::ostringstream oss;
            oss << " AFOLP (folp(0)=" << std::fixed << std::setprecision(3)
                << std::setw(6) << folp[0] << ")";
            s3 = oss.str();
        }

        title[ntitle] = s1 + s2 + s3;
        ntitle++;
    }
}

} // namespace feff::pot
