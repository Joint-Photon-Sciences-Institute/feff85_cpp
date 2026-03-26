// Many-pole self-energy grid calculation
// Converted from src/EXCH/mpse.f

#include "mpse.hpp"
#include "csigma.hpp"
#include "csigz.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include <feff/types.hpp>
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include "../common/logging.hpp"
#include "../common/file_io.hpp"

namespace feff::exch {

void mpse(const double edens[][252], int jintrs, FeffComplex e,
          int ipl, double mu)
{
    constexpr int NRPts = 12;

    double wpcorr[MxPole] = {};
    double gamma_arr[MxPole] = {};
    double ampfac[MxPole] = {};
    double delrHL[NRPts] = {};
    double deliHL[NRPts] = {};
    double rs1[NRPts] = {};

    FeffComplex zrnrm(0.0, 0.0);

    // Read exc.dat
    // Fortran: edens(jIntrs+1, 1) is the interstitial density for potential 1
    // In C++: edens[1][jintrs+1] (0-based potential index 1, 0-based radial index jintrs+1)
    double dens_int = edens[1][jintrs + 1];
    double rs_int = std::pow(3.0 / (4.0 * pi * dens_int), third);
    double wp_scale = std::sqrt(3.0 / (rs_int * rs_int * rs_int));

    int ipole = 0;
    {
        std::ifstream fin("exc.dat");
        feff::common::check_file_open(fin, "exc.dat", "mpse");

        std::string line;
        for (ipole = 0; ipole < MxPole; ++ipole) {
            // Skip comment lines (lines starting with #, *, c, C)
            feff::common::skip_comments(fin, "#*cC");
            if (!(fin >> wpcorr[ipole] >> gamma_arr[ipole] >> ampfac[ipole])) {
                break;
            }
            gamma_arr[ipole] /= hart;
            wpcorr[ipole] = (wpcorr[ipole] / hart) / wp_scale;
        }
    }
    // Sentinel
    if (ipole < MxPole) {
        wpcorr[ipole] = -1.0e30;
    }

    // Find the minimum and maximum densities across all potentials
    double max_dens = -1.0e30;
    double min_dens =  1.0e30;
    for (int ir = 0; ir < 251; ++ir) {
        for (int iph = 0; iph <= nphx; ++iph) {
            double d = edens[iph][ir];
            if (d > max_dens) max_dens = d;
            if (d < min_dens) min_dens = d;
        }
    }

    // Calculate RsMax(MinDens), RsMin(MaxDens)
    int rs_max = static_cast<int>(std::min(10.0, std::pow(3.0 / (4.0 * pi * min_dens), third)));
    int rs_min = static_cast<int>(std::max(0.001, std::pow(3.0 / (4.0 * pi * max_dens), third)));
    int rs_interst = static_cast<int>(std::pow(3.0 / (4.0 * pi * dens_int), third));

    // Calculate Sigma on a grid of NRPts points from RsMin to RsMax
    std::ofstream fout("mpse.bin");
    feff::common::check_file_open(fout, "mpse.bin", "mpse");

    double dr = static_cast<double>(rs_max - rs_min) / (NRPts - 1);

    for (int ir = 0; ir < NRPts; ++ir) {
        delrHL[ir] = 0.0;
        deliHL[ir] = 0.0;
        double rs_val = static_cast<double>(rs_min) + static_cast<double>(ir) * dr;
        rs1[ir] = rs_val;

        FeffComplex ztemp(0.0, 0.0);

        if (ipl > 1) {
            if (ipl == 2 || ir == NRPts - 1) {
                csigz(e, mu, rs_val, delrHL[ir], deliHL[ir], ztemp,
                      wpcorr, ampfac);
            } else if (ipl == 3) {
                delrHL[ir] = 0.0;
                deliHL[ir] = 0.0;
            } else {
                delrHL[ir] = delrHL[NRPts - 1];
            }
            if (ir == NRPts - 1) zrnrm = ztemp;
        } else {
            csigma(e, mu, rs1[ir], delrHL[ir], deliHL[ir], wpcorr, ampfac);
        }

        // Write Sigma(Rs, E) to mpse.bin
        fout << std::fixed << std::setprecision(10) << std::setw(20)
             << e.real() << std::setw(20) << rs_val
             << std::setw(20) << delrHL[ir] << std::setw(20) << deliHL[ir]
             << std::setw(20) << ztemp.real() << std::setw(20) << ztemp.imag()
             << "\n";
    }

    fout.close();
}

} // namespace feff::exch
