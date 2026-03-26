// Termination matrix M (energy-dependent, with polarization).
// Converted from GENFMT/mmtrxi.f

#include "mmtrxi.hpp"
#include <feff/constants.hpp>
#include <cmath>
#include <algorithm>

namespace feff::genfmt {

void mmtrxi(const FeffComplex rkk[][8], int lam1x,
            const BmatrixData& bmatd, int ie, int ileg, int ilegp,
            const int lind[8],
            const ClmzData& clmz, const LambdaData& lam,
            const NlmData& nlm, const double eta[],
            FmatrixData& fmat) {

    // Local gam and gamtl arrays
    FeffComplex gam_arr[ltot + 1][mtot + 1][ntot + 1];
    FeffComplex gamtl[ltot + 1][mtot + 1][ntot + 1];

    // Initialize
    for (int il = 0; il <= ltot; ++il)
        for (int im = 0; im <= mtot; ++im)
            for (int in = 0; in <= ntot; ++in) {
                gam_arr[il][im][in] = FeffComplex(0.0, 0.0);
                gamtl[il][im][in] = FeffComplex(0.0, 0.0);
            }

    // Set limits for orbital momentum from lind array
    int lmn = ltot;
    int lmx = 0;
    for (int k1 = 0; k1 < 8; ++k1) {
        if (lind[k1] > lmx) lmx = lind[k1];
        if (lind[k1] < lmn && lind[k1] >= 0) lmn = lind[k1];
    }
    int iln = lmn + 1;  // Fortran 1-based
    int ilx = lmx + 1;  // Fortran 1-based

    // Calculate factors gam and gamtl
    for (int il = iln; il <= ilx; ++il) {
        FeffComplex tltl_val = FeffComplex(static_cast<double>(2 * il - 1), 0.0);

        for (int l = 0; l < lam1x; ++l) {
            int m = lam.mlam[l];
            if (m < 0) continue;
            int im = m;  // 0-based
            if (im + 1 > il) continue;  // Fortran: if (im .gt. il)
            int in_val = lam.nlam[l];  // 0-based
            int imn_f = (in_val + 1) + m;  // Fortran 1-based imn = in + m

            int il0 = il - 1;  // 0-based

            double cam_sign = (m % 2 == 0) ? 1.0 : -1.0;
            FeffComplex cam = nlm.xnlm[il0][m] * cam_sign;

            if (imn_f <= il) {
                gam_arr[il0][m][in_val] = cam * clmz.clmi[il0][imn_f - 1][ileg];
            } else {
                gam_arr[il0][m][in_val] = FeffComplex(0.0, 0.0);
            }

            FeffComplex camt = tltl_val / nlm.xnlm[il0][m];
            gamtl[il0][m][in_val] = camt * clmz.clmi[il0][in_val][ilegp];
        }
    }

    // Compute fmati(lam1, lam2, ilegp) using bmati and rkk
    for (int lam1 = 0; lam1 < lam1x; ++lam1) {
        int m1 = lam.mlam[lam1];
        int in1 = lam.nlam[lam1];
        int iam1 = std::abs(m1);

        for (int lam2 = 0; lam2 < lam1x; ++lam2) {
            int m2 = lam.mlam[lam2];
            int in2 = lam.nlam[lam2];
            int iam2 = std::abs(m2);

            fmat.fmati[lam1][lam2][ilegp] = FeffComplex(0.0, 0.0);

            for (int k1 = 0; k1 < 8; ++k1) {
                for (int k2 = 0; k2 < 8; ++k2) {
                    int l1 = lind[k1] + 1;  // Fortran 1-based
                    int l2 = lind[k2] + 1;

                    if (l1 > 0 && l2 > 0 &&
                        iam1 + 1 <= l1 && iam2 + 1 <= l2) {
                        int l1_0 = l1 - 1;  // 0-based
                        int l2_0 = l2 - 1;

                        fmat.fmati[lam1][lam2][ilegp] +=
                            bmatd.bmati[m1 + mtot][k1][m2 + mtot][k2] *
                            rkk[ie][k1] * rkk[ie][k2] *
                            gam_arr[l1_0][iam1][in1] *
                            gamtl[l2_0][iam2][in2];
                    }
                }
            }

            // ileg is 0-based for clmz/dri but eta uses Fortran 1-based indexing
            fmat.fmati[lam1][lam2][ilegp] *=
                std::exp(-coni * eta[ileg + 1] * static_cast<double>(m1));
        }
    }
}

} // namespace feff::genfmt
