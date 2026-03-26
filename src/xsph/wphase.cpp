// Write diagnostic phase shift files.
// Converted from src/XSPH/wphase.f

#include "wphase.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace feff::xsph {

// Helper to index into ph array: ph[ie][ll+ltot][isp][iph]
static inline int ph_idx(int ie, int ll, int isp, int iph) {
    return ie + nex * ((ll + ltot) + (2 * ltot + 1) * (isp + nspx * iph));
}

// Helper to index into eref array: eref[ie][isp]
static inline int eref_idx(int ie, int isp) {
    return ie + nex * isp;
}

void wphase(int nph, const FeffComplex em[], const FeffComplex* eref,
            const int lmax[], int ne, const FeffComplex* ph,
            int ntitle, const char title[][80]) {

    for (int iph = 0; iph <= nph; iph++) {
        int linit = 0;
        if (linit >= lmax[iph] - 1) linit = lmax[iph] - 2;
        if (linit < 0) linit = 0;

        // Phase file
        std::ostringstream fname1;
        fname1 << "phase" << std::setfill('0') << std::setw(2) << iph << ".dat";
        std::ofstream fph(fname1.str());

        std::ostringstream fname2;
        fname2 << "phmin" << std::setfill('0') << std::setw(2) << iph << ".dat";
        std::ofstream fpm(fname2.str());

        for (int i = 0; i < ntitle; i++) {
            fph << "# " << title[i] << "\n";
            fpm << "# " << title[i] << "\n";
        }
        fph << "# " << iph << " " << lmax[iph] << " " << ne << "   unique pot,  lmax, ne\n";
        fpm << "# " << iph << " " << lmax[iph] << " " << ne << "   unique pot,  lmax, ne\n";
        fpm << "#      energy      re(eref)     re(p)    phase(" << linit
            << ")  phase(" << linit + 1 << ") phase(" << linit + 2 << ")\n";

        for (int ie = 0; ie < ne; ie++) {
            FeffComplex er = eref[eref_idx(ie, 0)];
            FeffComplex p2 = em[ie] - er;
            FeffComplex momentum = std::sqrt(2.0 * p2);

            fph << std::scientific << std::setprecision(6);
            fph << "#   ie        energy      re(eref)      im(eref)"
                << "         re(p)         im(p)\n";
            fph << std::setw(4) << ie
                << std::setw(14) << std::real(em[ie])
                << std::setw(14) << std::real(er)
                << std::setw(14) << std::imag(er)
                << std::setw(14) << std::real(momentum)
                << std::setw(14) << std::imag(momentum) << "\n";

            // Write phase shifts
            for (int ll = 0; ll <= lmax[iph]; ll++) {
                FeffComplex phval = ph[ph_idx(ie, ll, 0, iph)];
                fph << std::setw(14) << std::real(phval)
                    << std::setw(14) << std::imag(phval);
                if ((ll + 1) % 2 == 0 || ll == lmax[iph]) fph << "\n";
            }

            // Compact output
            fpm << std::scientific << std::setprecision(5)
                << std::setw(13) << std::real(em[ie])
                << std::setw(13) << std::real(er)
                << std::setw(13) << std::real(momentum);
            for (int ll = linit; ll <= linit + 2 && ll <= lmax[iph]; ll++) {
                fpm << std::setw(13) << std::real(ph[ph_idx(ie, ll, 0, iph)]);
            }
            fpm << "\n";
        }
    }
}

} // namespace feff::xsph
