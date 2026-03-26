#pragma once

// Termination matrix M (energy-independent part).
// Converted from GENFMT/mmtr.f
//
// Calculates the part of matrix M which does not depend on energy
// (see Rehr and Albers paper). For path expansion, always neglects
// spin-flip processes to simplify calculations.

#include "genfmt_data.hpp"

namespace feff::genfmt {

/// Calculate the energy-independent termination matrix B.
///
/// ipol: polarization flag
/// ispin: spin flag
/// le2: multipole moment selector
/// angks: angle between k-vector and spin-vector
/// ptz: polarization tensor (3x3, indexed [i+1][j+1] for i,j in {-1,0,1})
/// lind: orbital momentum indices for 8 transitions (output from bcoef)
/// kinit: initial kappa quantum number
/// ilinit: linit + 1
///
/// Output: bmat is filled.
void mmtr(BmatrixData& bmat, int ipol, int ispin, int le2, double angks,
          const FeffComplex ptz[3][3], int lind[8],
          const RotationMatrixData& rm, const double eta[],
          int nsc, int nleg, int kinit, int ilinit);

} // namespace feff::genfmt
