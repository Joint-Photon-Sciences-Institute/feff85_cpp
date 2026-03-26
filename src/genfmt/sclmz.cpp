// Spherical wave factors c_l^(m) z^m / m!
// Converted from GENFMT/sclmz.f
//
// Calculates energy dependent factors:
//   c(il,im) = c_l^(m) z^m / m! = c_lm  by recursion
//   c_{l+1,m} = c_{l-1,m} - (2l+1)*z*(c_{l,m} - c_{l,m-1})  (l != m)
//   c_{m,m} = (-z)^m * (2m)! / (2^m * m!)  with z = 1/(i*rho)

#include "sclmz.hpp"
#include <feff/constants.hpp>
#include <algorithm>

namespace feff::genfmt {

void sclmz(const FeffComplex rho[], int lmaxp1, int mmaxp1, int ileg,
           ClmzData& clmz) {

    FeffComplex cmm(1.0, 0.0);
    FeffComplex z = -coni / rho[ileg];

    // Fortran: clmi(1,1,ileg) = (1,0) -> C++: clmi[0][0][ileg]
    clmz.clmi[0][0][ileg] = FeffComplex(1.0, 0.0);
    // clmi(2,1,ileg) = clmi(1,1,ileg) - z
    clmz.clmi[1][0][ileg] = clmz.clmi[0][0][ileg] - z;

    int lmax = lmaxp1 - 1;

    // Recursion for m=0: c_{l+1,0} = c_{l-1,0} - z*(2l-1)*c_{l,0}
    // Fortran loop: do il = 2, lmax
    for (int il = 2; il <= lmax; ++il) {
        // Fortran: clmi(il+1,1,ileg) = clmi(il-1,1,ileg) - z*(2*il-1)*clmi(il,1,ileg)
        // C++ 0-based: [il][0][ileg] = [il-2][0][ileg] - z*(2*il-1)*[il-1][0][ileg]
        clmz.clmi[il][0][ileg] =
            clmz.clmi[il - 2][0][ileg] -
            z * static_cast<double>(2 * il - 1) * clmz.clmi[il - 1][0][ileg];
    }

    int mmxp1 = std::min(mmaxp1, lmaxp1);

    // Recursion for m >= 1
    // Fortran: do im = 2, mmxp1  (im is 1-based index, m = im-1)
    for (int im = 2; im <= mmxp1; ++im) {
        int m = im - 1;
        int imp1 = im + 1;  // Fortran imp1 = im+1; used as loop bound for il, NOT an array index

        cmm = -cmm * static_cast<double>(2 * m - 1) * z;

        // Fortran: clmi(im, im, ileg) = cmm
        // C++ 0-based: clmi[im-1][im-1][ileg]
        clmz.clmi[im - 1][im - 1][ileg] = cmm;

        // Fortran: clmi(imp1, im, ileg) = cmm * (2*m+1) * (1 - im*z)
        // C++ 0-based: clmi[im][im-1][ileg]
        clmz.clmi[im][im - 1][ileg] =
            cmm * static_cast<double>(2 * m + 1) * (1.0 - static_cast<double>(im) * z);

        // Fortran: do il = imp1, lmax
        for (int il = imp1; il <= lmax; ++il) {
            int l = il - 1;
            // Fortran: clmi(il+1,im,ileg) = clmi(l,im,ileg)
            //          - (2*l+1)*z*(clmi(il,im,ileg) + clmi(il,m,ileg))
            // C++ 0-based: [il][im-1][ileg] = [l-1][im-1][ileg]
            //          - (2*l+1)*z*([il-1][im-1][ileg] + [il-1][m-1][ileg])
            clmz.clmi[il][im - 1][ileg] =
                clmz.clmi[l - 1][im - 1][ileg] -
                static_cast<double>(2 * l + 1) * z *
                (clmz.clmi[il - 1][im - 1][ileg] + clmz.clmi[il - 1][m - 1][ileg]);
        }
    }
}

} // namespace feff::genfmt
