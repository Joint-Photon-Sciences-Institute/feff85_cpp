#pragma once

// Von Barth-Hedin exchange-correlation potential
// Converted from src/EXCH/vbh.f
// Reference: Von Barth, Hedin, J.Phys.C, 5, 1629, (1972), eq.6.2

namespace feff::exch {

// Von Barth-Hedin XC potential.
// Input:  rs    - density parameter
//         xmag  - 2 * fraction of given spin orientation
// Output: vxc   - XC potential for given spin orientation (Hartrees)
void vbh(double rs, double xmag, double& vxc);

// Auxiliary function for vbh.
// flarge(x) = (1 + x^3)*ln(1 + 1/x) + x/2 - x^2 - 1/3
double flarge(double x);

} // namespace feff::exch
