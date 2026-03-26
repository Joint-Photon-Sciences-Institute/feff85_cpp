#pragma once

// Termination matrix M (energy-dependent, with polarization).
// Converted from GENFMT/mmtrxi.f
//
// Calculates matrix M in the Rehr-Albers paper for the termination
// of the F-matrix product, incorporating reduced matrix elements rkk.

#include "genfmt_data.hpp"

namespace feff::genfmt {

/// Calculate the termination matrix M (fills fmat.fmati[...][...][ilegp]).
///
/// rkk: reduced matrix elements rkk[ie][kdif] (nex x 8)
/// lam1x: limit on lambda
/// bmat: energy-independent B matrix
/// ie: energy index (0-based)
/// ileg: leg index (0-based)
/// ilegp: leg' index (0-based), destination in fmati
/// lind: orbital momentum indices (8 transitions)
void mmtrxi(const FeffComplex rkk[][8], int lam1x,
            const BmatrixData& bmat, int ie, int ileg, int ilegp,
            const int lind[8],
            const ClmzData& clmz, const LambdaData& lam,
            const NlmData& nlm, const double eta[],
            FmatrixData& fmat);

} // namespace feff::genfmt
