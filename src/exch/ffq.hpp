#pragma once

// Analytic integral for HL imaginary self-energy
// Converted from src/EXCH/ffq.f

namespace feff::exch {

// Analytic integral used in imhl.
// Input:  q    - dimensionless momentum (normalized to Fermi momentum)
//         ef   - Fermi energy
//         xk   - momentum (inv Bohrs)
//         wp   - plasma frequency (units of Fermi energy)
//         alph - Hedin-Lundqvist parameter
// Output: return value
double ffq(double q, double ef, double xk, double wp, double alph);

} // namespace feff::exch
