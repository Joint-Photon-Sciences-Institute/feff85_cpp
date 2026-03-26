#pragma once

// Dirac-Hara exchange potential (energy-dependent exchange-correlation)
// Converted from src/EXCH/edp.f
// Reference: S.H. Chou, J.J. Rehr, E.A. Stern, E.R. Davidson (1986)

namespace feff::exch {

// Dirac-Hara exchange potential.
// Input:  rs - density parameter (a.u.)
//         xk - momentum (a.u.)
// Output: vr - Dirac potential (Hartrees)
void edp(double rs, double xk, double& vr);

} // namespace feff::exch
