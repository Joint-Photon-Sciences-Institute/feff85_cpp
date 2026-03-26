// Path importance factor calculation.
// Converted from GENFMT/import.f

#include "import_calc.hpp"
#include "../math/integration.hpp"
#include <feff/constants.hpp>
#include <cmath>
#include <complex>

namespace feff::genfmt {

void import_calc(int ne1, int nsp, int ik0, double reff, double deg,
                 const double ckmag[], const FeffComplex em[],
                 const FeffComplex eref2[][nspx], const FeffComplex cchi[],
                 double& xportx, double& crit) {

    double ffmag[nex];
    FeffComplex ck[nex];
    FeffComplex eref[nex];

    // Make importance factor: deg * integral(|chi| * d|p|)
    for (int ie = 0; ie < ne1; ++ie) {
        // Get ck(ie) set correctly
        ck[ie] = std::sqrt(2.0 * (em[ie] - eref2[ie][0]));
        if (nsp == 2) {
            eref[ie] = (eref2[ie][0] + eref2[ie][nspx - 1]) / 2.0;
            ck[ie] = std::sqrt(2.0 * (em[ie] - eref[ie]));
        }
        FeffComplex ckp = ck[ie];
        double xlam0 = ck[ie].imag() - ckp.imag();
        ffmag[ie] = std::abs(cchi[ie] * std::exp(2.0 * reff * xlam0));
    }

    // Integrate from edge (ik0) to ne1
    int nemax = ne1 - ik0;
    double xport = 0.0;
    math::trap(&ckmag[ik0], &ffmag[ik0], nemax, xport);
    xport = std::abs(deg * xport);

    if (xportx <= 0.0) xportx = xport;
    crit = 100.0 * xport / xportx;
}

} // namespace feff::genfmt
