#pragma once

// Broadened plasmon Hedin-Lundqvist self-energy
// Converted from src/EXCH/rhlbp.f
// Uses interpolation for both real and imaginary parts from bphl.dat

#include <string>

namespace feff::exch {

// Broadened plasmon HL self-energy.
// Input:  rs - density parameter
//         xk - momentum (a.u.)
// Output: erl - real part of self-energy
//         eim - imaginary part of self-energy
// Note: reads bphl.dat on first call (cached for subsequent calls).
void rhlbp(double rs, double xk, double& erl, double& eim);

} // namespace feff::exch
