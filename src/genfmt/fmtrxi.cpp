// F-matrix calculation for scattering amplitude.
// Converted from GENFMT/fmtrxi.f

#include "fmtrxi.hpp"
#include <feff/constants.hpp>
#include <cmath>
#include <algorithm>

namespace feff::genfmt {

void fmtrxi(int lam1x, int lam2x, int ie, int ileg, int ilegp,
            const ClmzData& clmz, const LambdaData& lam,
            const NlmData& nlm, const RotationMatrixData& rm,
            const PhaseData& pd, FmatrixData& fmat) {

    // Local arrays for gam and gamtl factors
    FeffComplex gam[ltot + 1][mtot + 1][ntot + 1];
    FeffComplex gamtl[ltot + 1][mtot + 1][ntot + 1];

    // Initialize
    for (int il = 0; il <= ltot; ++il)
        for (int im = 0; im <= mtot; ++im)
            for (int in = 0; in <= ntot; ++in) {
                gam[il][im][in] = FeffComplex(0.0, 0.0);
                gamtl[il][im][in] = FeffComplex(0.0, 0.0);
            }

    // Calculate factors gam and gamtl
    // Fortran: iln=1, ilx = lmax(ie, ipot(ilegp)) + 1
    // ileg and ilegp are 0-based for clmz/dri arrays, but ipot and eta
    // use the original Fortran 1-based convention, so we add 1.
    int ilegp_f = ilegp + 1;  // Fortran 1-based ilegp for ipot/eta access
    int ileg_f = ileg + 1;    // Fortran 1-based ileg for eta access
    int iln = 1;  // Fortran 1-based
    int ilx = pd.lmax[ie][pd.ipot[ilegp_f]] + 1;  // Fortran 1-based

    for (int il = iln; il <= ilx; ++il) {
        int ll = il - 1;

        // j-average t-matrix: tl = (exp(2i*ph(-ll)) - 1)/(2i)
        FeffComplex tl = (std::exp(2.0 * coni * ph_at(pd, ie, -ll, pd.ipot[ilegp_f])) - 1.0)
                         / (2.0 * coni);
        FeffComplex tltl = (std::exp(2.0 * coni * ph_at(pd, ie, ll, pd.ipot[ilegp_f])) - 1.0)
                           / (2.0 * coni);
        tltl = tl * static_cast<double>(ll + 1) + tltl * static_cast<double>(ll);

        int lam12x = std::max(lam1x, lam2x);

        for (int l = 0; l < lam12x; ++l) {
            int m = lam.mlam[l];
            if (m < 0) continue;
            int im = m;  // 0-based index (Fortran im = m+1, so our im = m)
            if (im + 1 > il) continue;  // Fortran: if (im .gt. il) goto 20
            int in = lam.nlam[l];  // 0-based (Fortran in = nlam(lam)+1, we use nlam[l])
            int imn = in + m;  // Fortran imn = in + m (both 1-based) = (nlam+1)+(m+1)-1

            // Adjust for Fortran 1-based indexing:
            // Fortran: gam(il, im, in) where im=m+1, in=nlam(lam)+1
            // C++ 0-based: gam[il-1][m][nlam[l]]
            int il0 = il - 1;  // 0-based il
            int in0 = lam.nlam[l];  // 0-based in (Fortran in = nlam(lam)+1)
            int imn_f = (in0 + 1) + m;  // Fortran imn = in + m (1-based)

            if (l < lam1x) {
                double cam_sign = (m % 2 == 0) ? 1.0 : -1.0;
                FeffComplex cam = nlm.xnlm[il0][m] * cam_sign;
                // Fortran: if (imn .le. il) gam(il,im,in) = cam * clmi(il,imn,ileg)
                // Fortran indices: il (1-based), imn (1-based), ileg (1-based)
                // C++: clmz.clmi[il-1][imn-1][ileg]
                if (imn_f <= il) {
                    gam[il0][m][in0] = cam * clmz.clmi[il0][imn_f - 1][ileg];
                } else {
                    gam[il0][m][in0] = FeffComplex(0.0, 0.0);
                }
            }

            if (l < lam2x) {
                FeffComplex camt = tltl / nlm.xnlm[il0][m];
                // Fortran: gamtl(il,im,in) = camt * clmi(il,in,ilegp)
                // C++: clmz.clmi[il-1][in_f-1][ilegp] where in_f = nlam(lam)+1
                gamtl[il0][m][in0] = camt * clmz.clmi[il0][in0][ilegp];
            }
        }
    }

    // Compute fmati(lam1, lam2, ilegp)
    for (int lam1 = 0; lam1 < lam1x; ++lam1) {
        int m1 = lam.mlam[lam1];
        int in1 = lam.nlam[lam1];  // 0-based (Fortran in1 = nlam(lam1)+1)
        int iam1 = std::abs(m1);
        int imn1_f = iam1 + 1 + in1;  // Fortran: imn1 = iam1 + in1 - 1 (with 1-based im,in)
        // Actually: Fortran iam1 = abs(m1)+1 (1-based), in1 = nlam(lam1)+1 (1-based)
        // imn1 = iam1 + in1 - 1 = abs(m1)+1 + nlam+1 - 1 = abs(m1) + nlam + 1
        // For 0-based: iam1_0 = abs(m1), in1_0 = nlam[lam1]

        for (int lam2 = 0; lam2 < lam2x; ++lam2) {
            int m2 = lam.mlam[lam2];
            int in2 = lam.nlam[lam2];
            int iam2 = std::abs(m2);

            FeffComplex cterm(0.0, 0.0);

            // Fortran: ilmin = max(iam1, iam2, imn1, in2, iln) -- all 1-based
            // In 1-based: iam1_f = abs(m1)+1, iam2_f = abs(m2)+1
            // in2_f = nlam[lam2]+1, imn1_f = abs(m1)+nlam[lam1]+1, iln=1
            int iam1_f = iam1 + 1;
            int iam2_f = iam2 + 1;
            int in2_f = in2 + 1;
            int imn1_fortran = iam1 + in1 + 1;
            int ilmin = std::max({iam1_f, iam2_f, imn1_fortran, in2_f, iln});

            for (int il = ilmin; il <= ilx; ++il) {
                // Fortran: skip terms with mu > l (NB il=l+1, so mu=il is mu>l)
                if (std::abs(m1) >= il || std::abs(m2) >= il) continue;

                // Rotation matrix: dri(il, m1d, m2d, ilegp) (Fortran 1-based)
                // m1d = m1 + mtot + 1 (Fortran), C++ index = m1 + mtot
                int il0 = il - 1;
                double dri_val = rm.dri[il0][m1 + mtot][m2 + mtot][ilegp];

                cterm += gam[il0][iam1][in1] * gamtl[il0][iam2][in2] * dri_val;
            }

            if (pd.eta[ileg_f] != 0.0) {
                cterm *= std::exp(-coni * pd.eta[ileg_f] * static_cast<double>(m1));
            }

            fmat.fmati[lam1][lam2][ilegp] = cterm;
        }
    }
}

} // namespace feff::genfmt
