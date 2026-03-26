#pragma once

// Polynomial root finding
// Converted from src/MATH/czeros.f (quadratic, cubic) and src/MATH/quartc.f

#include <feff/types.hpp>

namespace feff::math {

// Find zeros of quadratic polynomial: coef[0]*x^2 + coef[1]*x + coef[2] = 0
// Returns number of solutions in nsol (-1 if degenerate).
// Solutions stored in sol[0..nsol-1].
void cqdrtc(const FeffComplex coef[3], FeffComplex sol[2], int& nsol);

// Find zeros of cubic polynomial: coef[0]*x^3 + coef[1]*x^2 + coef[2]*x + coef[3] = 0
// Returns number of solutions in nsol (-1 if degenerate).
// Solutions stored in sol[0..nsol-1].
void ccubic(const FeffComplex coef[4], FeffComplex sol[4], int& nsol);

// Find roots of quartic polynomial: q[0]*x^4 + q[1]*x^2 + q[2]*x + q[3] = 0
// Note: no x^3 term. Input q[0..3] = {a, b, c, d}. Roots returned in q[0..3].
void quartic(FeffComplex q[4]);

} // namespace feff::math
