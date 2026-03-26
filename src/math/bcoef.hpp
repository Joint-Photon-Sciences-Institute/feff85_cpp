#pragma once

// Energy-independent B-matrix for absorption calculation
// Converted from src/MATH/bcoef.f
// Written by Alexei Ankudinov, March 2000

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <feff/types.hpp>

namespace feff::math {

// Calculate bmat: the energy-independent sum over polarization and
// angular momenta indices for absorption fine structure calculations.
//
// Input:
//   kinit  - kappa for initial orbital
//   ipol   - polarization type (0=average, 1=linear, 2=circular)
//   ptz    - polarization tensor (3x3, indexed -1:1,-1:1 → use [i+1][j+1])
//   le2    - multipole moment selector (0=E1 only, 1=M1, 2=E2)
//   ltrace - true for xsect.f (trace over ml)
//   ispin  - spin flag
//   angks  - angle between k-vector and spin-vector
//
// Output:
//   kind   - kappa values for 8 transitions
//   lind   - orbital momenta for 8 transitions
//   bmat   - energy independent matrix, dimensions:
//            bmat[ml+lx][ms][k-1][ml'+lx][ms'][k'-1]
//            where ml in [-lx,lx], ms in [0,1], k in [1,8]
void bcoef(int kinit, int ipol, const FeffComplex ptz[3][3], int le2,
           bool ltrace, int ispin, double angks,
           int kind[8], int lind[8],
           FeffComplex bmat[(2*lx+1)*2*8*(2*lx+1)*2*8]);

// Simplified interface using flat indexing
// bmat is indexed as bmat[ml1+lx][ms1][k1][ml2+lx][ms2][k2]
// Total size: (2*lx+1) * 2 * 8 * (2*lx+1) * 2 * 8
// Use bmat_index() to compute flat index.
inline int bmat_index(int ml1, int ms1, int k1, int ml2, int ms2, int k2) {
    // ml ranges [-lx, lx], ms ranges [0,1], k ranges [0,7]
    int i1 = (ml1 + lx) + (2 * lx + 1) * (ms1 + 2 * k1);
    int i2 = (ml2 + lx) + (2 * lx + 1) * (ms2 + 2 * k2);
    int dim1 = (2 * lx + 1) * 2 * 8;
    return i1 + dim1 * i2;
}

} // namespace feff::math
