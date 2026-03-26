// Overlap densities from neighbors.
// Converted from src/POT/ovrlp.f

#include "ovrlp.hpp"
#include "sumax.hpp"
#include "frnrm.hpp"
#include "../math/distance.hpp"
#include "../common/radial_grid.hpp"
#include <feff/constants.hpp>
#include <cmath>

namespace feff::pot {

void ovrlp(int iph, const int* iphat, const double* rat, const int* iatph,
           const int* novr, const int* iphovr, const int* nnovr, const double* rovr,
           const int* iz, int nat, const double* rho, double* dmag,
           const double* rhoval, const double* vcoul,
           double* edens, double* edenvl, double* vclap, double* rnrm)
{
    // Array access helpers (column-major, matching Fortran layout)
    // rho(i, iph) => rho[iph * 251 + i]  where iph is 0-based, i is 0-based
    // vcoul, rhoval, dmag: dimension (251, 0:nphx+1)  => stride 251, (nphx+2) columns
    // edens, edenvl, vclap: dimension (251, 0:nphx) => stride 251, (nphx+1) columns
    // iphovr(iovr, iph) => iphovr[iph * novrx + iovr]
    // nnovr(iovr, iph) => nnovr[iph * novrx + iovr]
    // rovr(iovr, iph) => rovr[iph * novrx + iovr]

    constexpr int s251 = 251;

    // Start with free atom values for current atom
    for (int i = 0; i < s251; ++i) {
        vclap[iph * s251 + i] = vcoul[iph * s251 + i];
        edens[iph * s251 + i] = rho[iph * s251 + i];
        edenvl[iph * s251 + i] = rhoval[iph * s251 + i];
    }

    if (novr[iph] > 0) {
        // Explicit overlap from OVERLAP card
        for (int iovr = 0; iovr < novr[iph]; ++iovr) {
            double rnn = rovr[iph * novrx + iovr];
            double ann = static_cast<double>(nnovr[iph * novrx + iovr]);
            int infr = iphovr[iph * novrx + iovr];
            sumax(rnn, ann, &vcoul[infr * s251], &vclap[iph * s251]);
            sumax(rnn, ann, &rho[infr * s251],   &edens[iph * s251]);
            sumax(rnn, ann, &rho[infr * s251],   &edenvl[iph * s251]);
        }
    } else {
        // Do overlapping from geometry with model atom iat
        int iat = iatph[iph];  // 1-based atom index from Fortran

        // Overlap with all atoms within r overlap max (rlapx)
        // 12 au = 6.35 ang
        double rlapx = 12.0;

        for (int inat = 0; inat < nat; ++inat) {
            // Don't overlap atom with itself (iat is 1-based, inat is 0-based)
            if (inat == iat - 1) continue;

            // If neighbor is too far away, don't overlap it
            double rnn = feff::math::dist(&rat[inat * 3], &rat[(iat - 1) * 3]);
            if (rnn > rlapx) continue;

            int infr = iphat[inat];
            sumax(rnn, one, &vcoul[infr * s251], &vclap[iph * s251]);
            sumax(rnn, one, &rho[infr * s251],   &edens[iph * s251]);
            sumax(rnn, one, &rho[infr * s251],   &edenvl[iph * s251]);
        }
    }

    // Set Norman radius
    frnrm(&edens[iph * s251], iz[iph], rnrm[iph]);

    // Remember ratio dmag/edens, not dmag itself
    for (int i = 0; i < s251; ++i) {
        if (edens[iph * s251 + i] > 0.0) {
            dmag[iph * s251 + i] = dmag[iph * s251 + i] / edens[iph * s251 + i];
        } else {
            dmag[iph * s251 + i] = 0.0;
        }
    }
}

} // namespace feff::pot
