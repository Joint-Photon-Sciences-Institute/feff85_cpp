#pragma once

// Wigner 3j coefficients and rotation matrices
// Converted from src/MATH/cwig3j.f and src/MATH/rotwig.f

namespace feff::math {

// Wigner 3j coefficient for integers (ient=1) or semi-integers (ient=2).
// Arguments j1,j2,j3,m1,m2 should be multiplied by ient.
double cwig3j(int j1, int j2, int j3, int m1, int m2, int ient);

// Wigner rotation matrix element using Wigner formula (Messiah eq.C.72).
// For integers (ient=1) or semi-integers (ient=2).
// Arguments jj,m1,m2 should be multiplied by ient.
double rotwig(double beta, int jj, int m1, int m2, int ient);

} // namespace feff::math
