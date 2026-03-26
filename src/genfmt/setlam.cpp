// Lambda parameter setup for genfmt.
// Converted from GENFMT/setlam.f

#include "setlam.hpp"
#include "../common/logging.hpp"
#include <feff/constants.hpp>
#include <cmath>
#include <cstdlib>
#include <sstream>

namespace feff::genfmt {

void setlam(int icalc, int ie, const double beta[], int nsc, int nleg,
            int ilinit, LambdaData& lam) {

    // One degree in radians
    constexpr double onedeg = 0.01745329252;

    int nmax = 0;
    int mmax = 0;
    int iord = 0;

    // Set iord, nmax and mmax based on icalc
    if (icalc < 0) {
        // Decode it and do what user wants
        int icode = -icalc;
        nmax = icode % 100;
        mmax = (icode % 10000) / 100;
        iord = icode / 10000 - 1;
    } else if (nsc == 1) {
        // Single scattering: exact
        mmax = ilinit;
        nmax = ilinit;
        iord = 2 * nmax + mmax;
    } else if (icalc < 10) {
        iord = icalc;
        mmax = iord;
        nmax = iord / 2;
    } else if (icalc == 10) {
        // Cute algorithm
        // Set mmax = L0 if straight line path, otherwise set mmax = 3
        mmax = ilinit;
        for (int ileg = 0; ileg < nleg; ++ileg) {
            int mag1 = static_cast<int>(std::abs(beta[ileg]));
            int mag2 = static_cast<int>(std::abs(mag1 - pi));
            // If beta is not 0 or pi, path is non-linear
            if (mag1 > onedeg && mag2 > onedeg) mmax = 3;
        }
        // Set nmax based on ie and l0
        // k <= 12 invA (ie=41)  nmax = L0
        // k >= 13 invA (ie=42)  nmax =  9
        nmax = ilinit;
        if (ie >= 42) nmax = 9;
        iord = 2 * nmax + mmax;
    } else {
        std::ostringstream ss;
        ss << " undefined icalc " << icalc;
        common::logger().wlog(ss.str());
        throw std::runtime_error("setlam: undefined icalc");
    }

    // Construct index lambda (lam), (mu, nu) = mlam(lam), nlam(lam)
    // Use mlam0/nlam0 for making indices, then sort
    int mlam0[lamtot];
    int nlam0[lamtot];

    int lam_count = 0;
    for (int in = 0; in <= nmax; ++in) {
        int n = in;
        for (int im = 0; im <= mmax; ++im) {
            int m = im;
            int jord = 2 * n + m;
            if (jord > iord) break;
            if (lam_count >= lamtot) {
                common::logger().wlog(" Lambda array filled, some order lost");
                goto done_fill;
            }
            mlam0[lam_count] = -m;
            nlam0[lam_count] = n;
            lam_count++;

            if (m == 0) continue;
            if (lam_count >= lamtot) {
                common::logger().wlog(" Lambda array filled, some order lost");
                goto done_fill;
            }
            mlam0[lam_count] = m;
            nlam0[lam_count] = n;
            lam_count++;
        }
    }
done_fill:
    lam.lamx = lam_count;

    if (lam.lamx > lamtot) {
        throw std::runtime_error("SETLAM lamx > lamtot");
    }

    // Sort mlam0 and nlam0 to use min possible laml0x
    // First pass: elements with n <= ilinit and |m| <= ilinit
    int sorted_count = 0;
    for (int lam0 = 0; lam0 < lam.lamx; ++lam0) {
        if (nlam0[lam0] <= ilinit &&
            std::abs(mlam0[lam0]) <= ilinit) {
            lam.nlam[sorted_count] = nlam0[lam0];
            lam.mlam[sorted_count] = mlam0[lam0];
            sorted_count++;
            nlam0[lam0] = -1;  // Mark as used
        }
    }
    lam.laml0x = sorted_count;

    // Second pass: remaining elements
    for (int lam0 = 0; lam0 < lam.lamx; ++lam0) {
        if (nlam0[lam0] >= 0) {
            lam.nlam[sorted_count] = nlam0[lam0];
            lam.mlam[sorted_count] = mlam0[lam0];
            sorted_count++;
        }
    }

    // Compute mmaxp1 and nmax from the sorted arrays
    lam.mmaxp1 = 0;
    lam.nmax = 0;
    for (int l = 0; l < lam.lamx; ++l) {
        if (lam.mlam[l] + 1 > lam.mmaxp1) lam.mmaxp1 = lam.mlam[l] + 1;
        if (lam.nlam[l] > lam.nmax) lam.nmax = lam.nlam[l];
    }

    if (lam.nmax > ntot || lam.mmaxp1 > mtot + 1) {
        std::ostringstream ss;
        ss << " mmaxp1, nmax, mtot, ntot "
           << lam.mmaxp1 << " " << lam.nmax << " " << mtot << " " << ntot;
        common::logger().wlog(ss.str());
        ss.str("");
        ss << " icalc " << icalc;
        common::logger().wlog(ss.str());
        throw std::runtime_error("setlam: nmax or mmaxp1 out of range");
    }
}

} // namespace feff::genfmt
