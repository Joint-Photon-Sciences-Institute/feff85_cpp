#pragma once

// Lambda parameter setup for genfmt.
// Converted from GENFMT/setlam.f

#include "genfmt_data.hpp"

namespace feff::genfmt {

/// Set lambda array based on calculation order and path type.
///
/// icalc: order of approximation
///   0  i0, ss exact
///   1  i1, ss exact
///   2  i2, ss exact
///  10  cute algorithm (energy- and geometry-dependent)
///  <0  explicit decode: icalc = -(nmax + 100*mmax + 10000*(iord+1))
///
/// ie: energy point index (used for cute algorithm)
/// nsc: number of scatterers
/// nleg: number of legs
/// ilinit: linit + 1
///
/// Outputs are stored in LambdaData.
void setlam(int icalc, int ie, const double beta[], int nsc, int nleg,
            int ilinit, LambdaData& lam);

} // namespace feff::genfmt
