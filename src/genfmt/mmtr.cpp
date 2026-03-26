// Termination matrix M (energy-independent part).
// Converted from GENFMT/mmtr.f

#include "mmtr.hpp"
#include "../math/bcoef.hpp"
#include <feff/constants.hpp>
#include <cmath>
#include <algorithm>

namespace feff::genfmt {

void mmtr(BmatrixData& bmatd, int ipol, int ispin, int le2, double angks,
          const FeffComplex ptz[3][3], int lind[8],
          const RotationMatrixData& rm, const double eta[],
          int nsc, int nleg, int kinit, int ilinit) {

    // Zero bmati
    for (int i = 0; i < 2 * mtot + 1; ++i)
        for (int j = 0; j < 8; ++j)
            for (int k = 0; k < 2 * mtot + 1; ++k)
                for (int l = 0; l < 8; ++l)
                    bmatd.bmati[i][j][k][l] = FeffComplex(0.0, 0.0);

    // Call bcoef to get bmat, kind, lind
    int kind[8];
    // bmat dimensions: (2*lx+1) * 2 * 8 * (2*lx+1) * 2 * 8
    constexpr int bmat_size = (2 * lx + 1) * 2 * 8 * (2 * lx + 1) * 2 * 8;
    FeffComplex bmat_flat[bmat_size];
    bool ltrace = false;

    math::bcoef(kinit, ipol, ptz, le2, ltrace, ispin, angks,
                kind, lind, bmat_flat);

    // Helper to access bmat_flat
    auto bmat = [&](int ml1, int ms1, int k1, int ml2, int ms2, int k2) -> FeffComplex {
        return bmat_flat[math::bmat_index(ml1, ms1, k1, ml2, ms2, k2)];
    };

    int is = 0;
    if (ispin == 1) is = nspx - 1;

    // ilinit = initial orb. momentum + 1
    int lxx = std::min(mtot, ilinit);

    // Set indices for bmati
    for (int mu1 = -lxx; mu1 <= lxx; ++mu1) {
        int mu1d = mu1 + mtot;  // 0-based index into dri array (Fortran: mu1+mtot+1)
        for (int mu2 = -lxx; mu2 <= lxx; ++mu2) {
            int mu2d = mu2 + mtot;

            if (ipol != 0) {
                for (int k1 = 0; k1 < 8; ++k1) {
                    for (int k2 = 0; k2 < 8; ++k2) {
                        int l1 = lind[k1];  // Already 0-based orbital momentum
                        int l2 = lind[k2];

                        for (int m1 = -lind[k1]; m1 <= lind[k1]; ++m1) {
                            for (int m2 = -lind[k2]; m2 <= lind[k2]; ++m2) {
                                int m1d = m1 + mtot;
                                int m2d = m2 + mtot;

                                // bmati(mu1, k1, mu2, k2) += bmat(m1,is,k1, m2,is,k2)
                                //   * exp(-i*(eta(nsc+2)*m2 + eta(0)*m1))
                                //   * dri(l1+1, mu1d, m1d, nsc+2)
                                //   * dri(l2+1, m2d, mu2d, nleg)
                                // Fortran: nsc+2 is the extra rotation index for polarization
                                // dri indices: l1+1 in Fortran 1-based = l1 in C++ 0-based
                                bmatd.bmati[mu1 + mtot][k1][mu2 + mtot][k2] +=
                                    bmat(m1, is, k1, m2, is, k2) *
                                    std::exp(-coni * (eta[nsc + 2] * static_cast<double>(m2) +
                                                      eta[0] * static_cast<double>(m1))) *
                                    rm.dri[l1][mu1d][m1d][nsc + 1] *
                                    rm.dri[l2][m2d][mu2d][nleg - 1];
                            }
                        }
                    }
                }
            } else {
                // ipol=0: bmat is diagonal in k1,k2 and LS L'S'
                // Two rotation matrices can be combined to 1
                for (int k1 = 0; k1 < 8; ++k1) {
                    int l1 = lind[k1];
                    if (l1 >= 0) {
                        int m1 = 0;
                        int m1d = m1 + mtot;

                        // dri(l1+1, mu1d, mu2d, nsc+1) in Fortran
                        // C++ 0-based: dri[l1][mu1d][mu2d][nsc]
                        bmatd.bmati[mu1 + mtot][k1][mu2 + mtot][k1] +=
                            bmat(m1, is, k1, m1, is, k1) *
                            rm.dri[l1][mu1d][mu2d][nsc];
                    }
                }
            }
        }
    }
}

} // namespace feff::genfmt
