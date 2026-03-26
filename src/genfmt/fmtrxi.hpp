#pragma once

// F-matrix calculation for scattering amplitude.
// Converted from GENFMT/fmtrxi.f
//
// f(lam,lam') = sum_l tl * gam(l,m,n) * dri(l,m,m',ileg) * gamt(l,m',n')
//             * exp(-i*m*eta)

#include "genfmt_data.hpp"

namespace feff::genfmt {

/// Calculate scattering amplitude matrix f(lam,lam') for a given leg pair.
///
/// lam1x, lam2x: limits on lambda and lambda'
/// ie: energy index (0-based)
/// ileg: leg index (0-based)
/// ilegp: leg' index (0-based), where the scatterer potential is used
///
/// Uses: clmz, lambda arrays, rotation matrices, phase shifts, eta
/// Output: fmat.fmati[...][...][ilegp] is set for the current energy point.
void fmtrxi(int lam1x, int lam2x, int ie, int ileg, int ilegp,
            const ClmzData& clmz, const LambdaData& lam,
            const NlmData& nlm, const RotationMatrixData& rm,
            const PhaseData& pd, FmatrixData& fmat);

} // namespace feff::genfmt
