#pragma once

// Associated Legendre polynomials
// Converted from src/MATH/cpl0.f

namespace feff::math {

// Calculate associated Legendre polynomials P_l0(x) by recursion.
// pl0 must have size >= lmaxp1. Output in pl0[0..lmaxp1-1].
void cpl0(double x, double pl0[], int lmaxp1);

} // namespace feff::math
