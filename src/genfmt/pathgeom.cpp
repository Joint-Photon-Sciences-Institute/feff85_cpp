// Path geometry computation (Euler angles and leg lengths).
// Converted from GENFMT/pathgeom.f

#include "pathgeom.hpp"
#include "trig_utils.hpp"
#include <feff/constants.hpp>
#include <cmath>

namespace feff::genfmt {

void pathgeom(int nleg, int& nsc, int ipol,
              double rat[][legtot + 2], int ipot[],
              double ri[], double beta[], double eta[]) {

    double alpha[legtot + 2];  // alpha(0:legtot+1)
    double gamma[legtot + 1];  // gamma(1:legtot+1) -> 0-based: gamma[0..legtot]

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

    // Fortran: do j = 1, nangle (1-based)
    for (int j = 1; j <= nangle; ++j) {
        int ifix = 0;
        int i_atom, ip1, im1;

        if (j == nsc + 1) {
            // j+1 'z' atom, j central atom, j-1 last path atom
            i_atom = 0;
            ip1 = 1;
            if (ipol > 0) ip1 = nleg + 1;
            im1 = nsc;
        } else if (j == nsc + 2) {
            // j central atom, j+1 first path atom, j-1 'z' atom
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

        // Handle special case
        if (ifix == 1) {
            x = 0.0;
            y = 0.0;
            z = 1.0;
            trig(x, y, z, ct, st, cp, sp);
            ifix = 0;
        }

        // cppp = cos(phi' - phi), sppp = sin(phi' - phi)
        double cppp = cp * cpp + sp * spp;
        double sppp = spp * cp - cpp * sp;
        double phi = std::atan2(sp, cp);
        double phip = std::atan2(spp, cpp);

        // alph = exp(i*alpha), beta = cos(beta), gamm = exp(i*gamma)
        FeffComplex alph = -(st * ctp - ct * stp * cppp - coni * stp * sppp);
        double beta_cos = ct * ctp + st * stp * cppp;
        // Watch out for roundoff
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

        alpha[j] = alpha_j;
        gamma[j] = gamma_j;

        if (j <= nleg) {
            // rat layout is [xyz][atom_index], so xyz are not contiguous per atom.
            double dx = rat[0][i_atom] - rat[0][im1];
            double dy = rat[1][i_atom] - rat[1][im1];
            double dz = rat[2][i_atom] - rat[2][im1];
            ri[j] = std::sqrt(dx * dx + dy * dy + dz * dz);
        }
    }

    // Make eta(i) = alpha(i-1) + gamma(i)
    // We need alpha(nangle) = alpha(0)
    alpha[0] = alpha[nangle];
    for (int j = 1; j <= nleg; ++j) {
        eta[j] = alpha[j - 1] + gamma[j];
    }
    if (ipol > 0) {
        eta[0] = gamma[nleg + 1];
        eta[nleg + 1] = alpha[nleg];
    }
}

} // namespace feff::genfmt
