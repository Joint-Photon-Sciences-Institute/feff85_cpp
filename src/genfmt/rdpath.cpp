// Read path data from paths.dat file.
// Converted from GENFMT/rdpath.f

#include "rdpath.hpp"
#include "trig_utils.hpp"
#include <feff/constants.hpp>
#include "../common/logging.hpp"
#include <cmath>
#include <sstream>
#include <string>

namespace feff::genfmt {

void rdpath(std::istream& in, bool& done, int ipol,
            std::string potlbl[], double rat[][legtot + 2],
            double ri[], double beta[], double eta[],
            double& deg, int ipot[], int& nsc, int& nleg,
            int npot, int& ipath) {

    std::string line;

    // Read: ipath, nleg, deg
    if (!std::getline(in, line)) {
        done = true;
        return;
    }

    {
        std::istringstream ss(line);
        if (!(ss >> ipath >> nleg >> deg)) {
            done = true;
            return;
        }
    }

    if (nleg > legtot) {
        std::ostringstream ss;
        ss << " nleg .gt. legtot, nleg, legtot " << nleg << " " << legtot;
        common::logger().wlog(ss.str());
        common::logger().wlog(" ERROR");
        done = true;
        return;
    }

    // Skip label line (x y z ipot rleg beta eta)
    std::getline(in, line);

    // Read atom coordinates and potential indices
    // Fortran: do ileg = 1, nleg (1-based)
    for (int ileg = 1; ileg <= nleg; ++ileg) {
        if (!std::getline(in, line)) {
            common::logger().wlog(" Unexpected end of file");
            throw std::runtime_error("rdpath: unexpected end of file");
        }
        std::istringstream ss(line);
        double rx, ry, rz;
        int ip;
        std::string plbl;
        ss >> rx >> ry >> rz >> ip >> plbl;

        // Convert from Angstrom to code units (Bohr)
        rat[0][ileg] = rx / bohr;
        rat[1][ileg] = ry / bohr;
        rat[2][ileg] = rz / bohr;
        ipot[ileg] = ip;
        if (!plbl.empty()) {
            potlbl[ip] = plbl;
        }

        if (ipot[ileg] > npot) {
            std::ostringstream es;
            es << " ipot(ileg) too big, ipot, ileg, npot "
               << ipot[ileg] << " " << ileg << " " << npot;
            common::logger().wlog(es.str());
            common::logger().wlog(" ERROR");
            done = true;
            return;
        }
    }

    nsc = nleg - 1;

    // We need the 'z' atom for polarization case
    if (ipol > 0) {
        rat[0][nleg + 1] = rat[0][nleg];
        rat[1][nleg + 1] = rat[1][nleg];
        rat[2][nleg + 1] = rat[2][nleg] + 1.0;
    }

    // Add rat(0) and ipot(0)
    for (int j = 0; j < 3; ++j) {
        rat[j][0] = rat[j][nleg];
    }
    ipot[0] = ipot[nleg];

    int nangle = nleg;
    if (ipol > 0) nangle = nleg + 1;

    double alpha_arr[legtot + 2];  // alpha(0:legtot+1)
    double gamma_arr[legtot + 2];  // gamma(1:legtot+1) -> stored at same indices

    for (int j = 1; j <= nangle; ++j) {
        int ifix = 0;
        int i_atom, ip1, im1;

        if (j == nsc + 1) {
            i_atom = 0;
            ip1 = 1;
            if (ipol > 0) ip1 = nleg + 1;
            im1 = nsc;
        } else if (j == nsc + 2) {
            i_atom = 0;
            ip1 = 1;
            im1 = nleg + 1;
            ifix = 1;
        } else {
            i_atom = j;
            ip1 = j + 1;
            im1 = j - 1;
        }

        double x, y, z;
        double ctp, stp, cpp, spp;
        double ct, st, cp, sp;

        x = rat[0][ip1] - rat[0][i_atom];
        y = rat[1][ip1] - rat[1][i_atom];
        z = rat[2][ip1] - rat[2][i_atom];
        trig(x, y, z, ctp, stp, cpp, spp);

        x = rat[0][i_atom] - rat[0][im1];
        y = rat[1][i_atom] - rat[1][im1];
        z = rat[2][i_atom] - rat[2][im1];
        trig(x, y, z, ct, st, cp, sp);

        if (ifix == 1) {
            x = 0.0;
            y = 0.0;
            z = 1.0;
            trig(x, y, z, ct, st, cp, sp);
            ifix = 0;
        }

        double cppp = cp * cpp + sp * spp;
        double sppp = spp * cp - cpp * sp;
        double phi = std::atan2(sp, cp);
        double phip = std::atan2(spp, cpp);

        FeffComplex alph = -(st * ctp - ct * stp * cppp - coni * stp * sppp);
        double beta_cos = ct * ctp + st * stp * cppp;
        if (beta_cos < -1.0) beta_cos = -1.0;
        if (beta_cos >  1.0) beta_cos =  1.0;
        FeffComplex gamm = -(st * ctp * cppp - ct * stp + coni * st * sppp);

        double alpha_j, gamma_j;
        arg(alph, phip - phi, alpha_j);
        beta[j] = std::acos(beta_cos);
        arg(gamm, phi - phi, gamma_j);

        // Convert from rotation of FRAME to rotation of VECTORS
        double dumm = alpha_j;
        alpha_j = pi - gamma_j;
        gamma_j = pi - dumm;

        alpha_arr[j] = alpha_j;
        gamma_arr[j] = gamma_j;

        if (j <= nleg) {
            double dx = rat[0][i_atom] - rat[0][im1];
            double dy = rat[1][i_atom] - rat[1][im1];
            double dz = rat[2][i_atom] - rat[2][im1];
            ri[j] = std::sqrt(dx * dx + dy * dy + dz * dz);
        }
    }

    // Make eta(i) = alpha(i-1) + gamma(i)
    alpha_arr[0] = alpha_arr[nangle];
    for (int j = 1; j <= nleg; ++j) {
        eta[j] = alpha_arr[j - 1] + gamma_arr[j];
    }
    if (ipol > 0) {
        eta[0] = gamma_arr[nleg + 1];
        eta[nleg + 1] = alpha_arr[nleg];
    }

    done = false;
}

} // namespace feff::genfmt
