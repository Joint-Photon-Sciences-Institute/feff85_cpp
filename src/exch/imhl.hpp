#pragma once

// Hedin-Lundqvist imaginary self-energy (analytic)
// Converted from src/EXCH/imhl.f
// Written by J. Mustre (March 1988), modified by J. Rehr (Oct 1991)

namespace feff::exch {

// HL imaginary self-energy.
// Input:  rs - density parameter
//         xk - momentum (a.u.)
// Output: eim   - imaginary self-energy
//         icusp - 1 if at cusp position, 0 otherwise
void imhl(double rs, double xk, double& eim, int& icusp);

} // namespace feff::exch
