#pragma once

// Quinn electron-hole pair losses
// Converted from src/EXCH/quinn.f
// Reference: J.J. Quinn, Phys. Rev. 126, 1453 (1962), eq. (7)

namespace feff::exch {

// Quinn approximation for e-h pair losses below plasmon turn on.
// Input:  x  - p/pf (momentum ratio)
//         rs - density parameter
//         wp - plasma frequency (units of Fermi energy)
//         ef - Fermi energy
// Output: ei - imaginary self energy
void quinn(double x, double rs, double wp, double ef, double& ei);

} // namespace feff::exch
