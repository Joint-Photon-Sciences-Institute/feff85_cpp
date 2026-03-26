#pragma once

// Spherical Bessel function calculations
// Converted from src/MATH/besjn.f, besjh.f, bjnser.f, exjlnl.f

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <feff/types.hpp>
#include <vector>

namespace feff::math {

// Calculate spherical Bessel functions jl and nl for l = 0 to ltot+1
// using series expansion, recursion, or asymptotic forms depending on |x|.
// jl and nl use Abramowitz conventions (nl = yl).
// Arrays must have size >= ltot+2.
void besjn(FeffComplex x, FeffComplex jl[], FeffComplex nl[]);

// Calculate spherical Bessel functions jl and hl for l = 0 to lbmax.
// hl = hl^+ (Messiah convention) for Im(x) >= 0
// hl = hl^- (Messiah convention) for Im(x) < 0
void besjh(FeffComplex x, int lbmax, FeffComplex jl[], FeffComplex hl[]);

// Series expansion for spherical Bessel functions jl and nl of order l.
// ifl: 0 = return both jl and nl, 1 = jl only, 2 = nl only
void bjnser(FeffComplex x, int l, FeffComplex& jl, FeffComplex& nl, int ifl);

// Exact analytic spherical Bessel functions jl and nl for l = 0 to 6.
// Falls back to series expansion for |z| < 0.3.
void exjlnl(FeffComplex z, int l, FeffComplex& jl, FeffComplex& nl);

} // namespace feff::math
