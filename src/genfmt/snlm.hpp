#pragma once

// Legendre normalization factors.
// Converted from GENFMT/snlm.f
//
// xnlm(l,m) = sqrt((2l+1)(l-m)!/(l+m)!)

#include <feff/dimensions.hpp>

namespace feff::genfmt {

/// Compute Legendre normalization factors xnlm(il,im) for il=0..lmaxp1-1, im=0..mmaxp1-1.
/// xnlm is indexed as xnlm[il][im] with dimensions (ltot+1, mtot+1).
void snlm(int lmaxp1, int mmaxp1, double xnlm[ltot + 1][mtot + 1]);

/// Compute scaled factorials: flg[i] = i! * afac^i.
/// Returns afac (scaling factor = 1/64).
double factst(double flg[211]);

} // namespace feff::genfmt
