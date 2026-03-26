#pragma once

// Convolution of cross-section with energy loss function
// Converted from src/MATH/conv.f

#include <feff/types.hpp>

namespace feff::math {

// Convolute xsec(omega) with Lorentzian loss function.
// omega: energy grid, xsec: cross-section (modified in place), ne1: number of points.
// vicorr: Lorentzian width parameter.
void conv(const double omega[], FeffComplex xsec[], int ne1, double vicorr);

} // namespace feff::math
