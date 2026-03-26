#pragma once
// Angular coefficient matrix for spin/orbital moment calculations.
// Converted from src/XSPH/acoef.f

#include <feff/dimensions.hpp>

namespace feff::xsph {

/// Compute angular coefficient matrix for density calculation.
///
/// @param ispin  Spin type: 0=occupation numbers, +-1=Sz/Lz/Tz, +-2=Sz/Nl/Nj
/// @param amat   Output matrix: amat[ml+lx][i1][i2][iop][lpp]
///               Dimensions: (2*lx+1) x 2 x 2 x 3 x (lx+1)
///               Stored as flat array, size (2*lx+1)*2*2*3*(lx+1)
void acoef(int ispin, float amat[]);

/// Helper: get kappa and j from index i and orbital momentum lpp.
/// @param i    Index (1 or 2)
/// @param lpp  Orbital angular momentum
/// @param jj   Output: total angular momentum - 1/2
/// @param k    Output: kappa quantum number (0 if invalid)
void kfromi(int i, int lpp, int& jj, int& k);

/// Index into amat array: amat[ml+lx][i1][i2][iop][lpp]
/// ml in [-lx,lx], i1 in [0,1], i2 in [0,1], iop in [0,2], lpp in [0,lx]
inline int amat_index(int ml, int i1, int i2, int iop, int lpp) {
    return (ml + lx) + (2 * lx + 1) * (i1 + 2 * (i2 + 2 * (iop + 3 * lpp)));
}

} // namespace feff::xsph
