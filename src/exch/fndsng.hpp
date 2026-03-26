#pragma once

// Singularity finder for Hedin-Lundqvist self-energy integrands
// Converted from src/EXCH/fndsng.f
//
// Finds singularities in the integrands of eq. 13 in:
//   "Single-particle Spectrum of the Degenerate Electron Gas
//    II. Numerical Results for Electrons Coupled to Plasmons"
//   Phys. kondens. Materie, Bd. 6 (1967)

#include <feff/types.hpp>

namespace feff::exch {

/// Find singularities in HL self-energy integrands between limit1 and limit2.
///
/// Solves for singularities of one of the three integrands, then checks
/// whether each singularity lies within the integration limits.
///
/// @param limit1  Lower limit of integration
/// @param limit2  Upper limit of integration
/// @param nsing   [out] Number of singularities found
/// @param xsing   [out] Array of singularities (max 20)
/// @param dppar   Double precision parameters:
///                dppar[0] = Wp/EFermi, dppar[1] = Gamma/EFermi,
///                dppar[2] = Energy/EFermi, dppar[3] = gap energy
/// @param cpar    Complex parameters:
///                cpar[0] = ck/kFermi, cpar[1] = Energy/EFermi + i*Gamma/EFermi
/// @param ifcn    Function index (1: solve eqs 1+2, 2: solve eq 1 only)
void fndsng(FeffComplex limit1, FeffComplex limit2,
            int& nsing, FeffComplex xsing[20],
            const double dppar[10], const FeffComplex cpar[10], int ifcn);

} // namespace feff::exch
