// Diagnostic potential output.
// Converted from src/POT/wpot.f

#include "wpot.hpp"
#include "../common/radial_grid.hpp"
#include "../common/logging.hpp"
#include "../common/file_io.hpp"
#include <feff/constants.hpp>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace feff::pot {

void wpot(int nph, const double* edens, const int* imt, const int* inrm,
          const double* rho, const double* vclap, const double* vcoul,
          const double* vtot, int ntitle, const std::string* title)
{
    using feff::common::rr;
    constexpr int s251 = 251;

    // Units: potentials in Hartrees (v*27.2 -> eV)
    //        density in #/(bohr)^3 (rho*e/(.529)^3 -> e/(Ang)^3)

    for (int iph = 0; iph <= nph; ++iph) {
        // Prepare filename: potXX.dat
        std::ostringstream fname;
        fname << "pot" << std::setw(2) << std::setfill('0') << iph << ".dat";

        std::ofstream fout(fname.str());
        feff::common::check_file_open(fout, fname.str(), "wpot");

        // Write header
        for (int it = 0; it < ntitle; ++it) {
            fout << " " << title[it] << "\n";
        }

        fout << " " << std::setw(3) << iph
             << std::setw(4) << imt[iph]
             << std::setw(4) << inrm[iph]
             << "  Unique potential, I_mt, I_norman."
             << "    Following data in atomic units.\n";

        fout << " iph " << iph << "\n";
        fout << "   i      r         vcoul        rho"
             << "     ovrlp vcoul  ovrlp vtot  ovrlp rho\n";

        // Write data up to r <= 38 (gives ~249 points)
        for (int i = 0; i < s251; ++i) {
            double r = rr(i + 1);  // 1-based Fortran grid
            if (r > 38.0) break;

            fout << " " << std::setw(4) << (i + 1)
                 << std::scientific << std::setprecision(4)
                 << std::setw(12) << r
                 << std::setw(12) << vcoul[iph * s251 + i]
                 << std::setw(12) << rho[iph * s251 + i] / (4.0 * pi)
                 << std::setw(12) << vclap[iph * s251 + i]
                 << std::setw(12) << vtot[iph * s251 + i]
                 << std::setw(12) << edens[iph * s251 + i] / (4.0 * pi)
                 << "\n";
        }
        fout.close();
    }
}

} // namespace feff::pot
