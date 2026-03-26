// Spherical overlap summation (Louck's method).
// Converted from src/POT/sumax.f

#include "sumax.hpp"
#include "../common/radial_grid.hpp"
#include <cmath>

namespace feff::pot {

void sumax(double rn, double ann, const double* aa2, double* aasum)
{
    // aa2 and aasum are 0-based (size 250), matching Fortran nptx=250 (1:250)
    constexpr int nptx = 250;
    double stor[nptx];

    using feff::common::ii;
    using feff::common::xx;

    // Wigner-Seitz radius is set to 15 in ATOM
    double rws = 15.0;
    int jjchi = ii(rws);   // 1-based index
    int jtop = ii(rn);     // 1-based index

    double topx = xx(jjchi);

    // xbl must persist across iterations, matching Fortran where it's initialized
    // once before the loop (xbl = 0) and then modified inside conditionally.
    double xbl = 0.0;

    // Fortran loop: do i = 1, jtop  (1-based)
    for (int i = 1; i <= jtop; ++i) {
        double x = xx(i);
        double xint = 0.0;
        double et = std::exp(x);
        double blx = std::log(rn - et);

        if (blx >= topx) {
            stor[i - 1] = 0.5 * xint * ann / (rn * et);
            continue;
        }

        int jbl = static_cast<int>(2.0 + 20.0 * (blx + 8.8));
        if (jbl < 1) jbl = 1;

        if (jbl >= 2) {
            // Use linear interp to make end cap near center of neighbor
            double xjbl = static_cast<double>(jbl);
            xbl = 0.05 * (xjbl - 1.0) - 8.8;
            double g = xbl - blx;
            // aa2 is 0-based: aa2[jbl-1] = Fortran aa2(jbl)
            xint = xint + 0.5 * g * (aa2[jbl - 1] * (2.0 - 20.0 * g) * std::exp(2.0 * xbl)
                   + 20.0 * g * aa2[jbl - 2] * std::exp(2.0 * (xbl - 0.05)));
        }

        double tlx = std::log(rn + et);
        int jtl;
        if (tlx >= topx) {
            jtl = jjchi;
        } else {
            jtl = static_cast<int>(1.0 + 20.0 * (tlx + 8.8));
            if (jtl < jbl) {
                // Handle peculiar special case at center of atom 1
                double fzn = aa2[jtl - 1] * std::exp(2.0 * (xbl - 0.05));
                double fz3 = aa2[jbl - 1] * std::exp(2.0 * xbl);
                double fz2 = fzn + 20.0 * (fz3 - fzn) * (tlx - xbl + 0.05);
                double fz1 = fzn + 20.0 * (fz3 - fzn) * (blx - xbl + 0.05);
                xint = 0.5 * (fz1 + fz2) * (tlx - blx);
                stor[i - 1] = 0.5 * xint * ann / (rn * et);
                continue;
            }
            double xjtl = static_cast<double>(jtl);
            double xtl = 0.05 * (xjtl - 1.0) - 8.8;
            double c = tlx - xtl;
            xint = xint + 0.5 * c * (aa2[jtl - 1] * (2.0 - 20.0 * c)
                   * std::exp(2.0 * xtl) + aa2[jtl] * 20.0 * c
                   * std::exp(2.0 * (xtl + 0.05)));
        }

        if (jtl > jbl) {
            // Trapezoidal sum in the interior
            // Must match Fortran: xbl is only incremented when jbl < jtl
            // after jbl++, so xbl is NOT incremented on the last iteration.
            // This matters because xbl persists across outer loop iterations.
            for (;;) {
                xint = xint + 0.5 * (aa2[jbl - 1] * std::exp(2.0 * xbl) +
                       aa2[jbl] * std::exp(2.0 * (xbl + 0.05))) * 0.05;
                jbl = jbl + 1;
                if (jbl < jtl) {
                    xbl = xbl + 0.05;
                } else {
                    break;
                }
            }
        }

        stor[i - 1] = 0.5 * xint * ann / (rn * et);
    }

    // Accumulate into aasum
    for (int i = 0; i < jtop; ++i) {
        aasum[i] = aasum[i] + stor[i];
    }
}

} // namespace feff::pot
