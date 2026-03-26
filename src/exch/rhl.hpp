#pragma once

// Hedin-Lundqvist self-energy (interpolated real + analytic imaginary)
// Converted from src/EXCH/rhl.f
// Written by Jose Mustre

namespace feff::exch {

// HL self-energy using interpolation for real part, analytic imaginary part.
// Input:  rs - density parameter
//         xk - momentum (a.u.)
// Output: erl - real part of self-energy
//         eim - imaginary part of self-energy
void rhl(double rs, double xk, double& erl, double& eim);

} // namespace feff::exch
