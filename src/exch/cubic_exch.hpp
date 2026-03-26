#pragma once

// Cubic dispersion solver for EXCH module
// Converted from src/EXCH/cubic.f
// Note: different from math/ccubic

namespace feff::exch {

// Finds roots of: 4*xk0 * q^3 + (alph - 4*xk0^2) * q^2 + wp^2 = 0
// See Abramowitz and Stegun pg 17 for formulae.
// Input:  xk0   - normalized momentum
//         wp    - plasma frequency
//         alph  - HL parameter
// Output: rad    - discriminant (q^3 + r^2)
//         qplus  - positive root
//         qminus - negative root
void cubic_exch(double xk0, double wp, double alph,
                double& rad, double& qplus, double& qminus);

} // namespace feff::exch
