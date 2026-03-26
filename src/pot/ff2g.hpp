#pragma once
// Integrate l-resolved DOS from FMS Green's function.
// Converted from src/POT/ff2g.f
//
// Combines the FMS chi (from gtr) with the embedded atom quantities
// to produce l-resolved DOS (xrhoce) and r-dependent DOS (yrhoce).

#include <feff/dimensions.hpp>
#include <feff/types.hpp>
#include <complex>

namespace feff::pot {

/// Integrate l-resolved and r-dependent DOS from FMS results.
///
/// @param gtr       FMS Green's function chi [lx+1] (single precision complex)
/// @param iph       Potential index
/// @param ie        Energy index (1-based)
/// @param ilast     Last radial grid index
/// @param xrhoce    Central atom l-DOS [(lx+1)][nphx+1] (flat, accumulated)
/// @param xrhole    Embedded atom l-DOS [lx+1] (input)
/// @param xrhocp    Previous-step central l-DOS [(lx+1)][nphx+1] (flat, accumulated)
/// @param ee        Current complex energy
/// @param ep        Previous complex energy
/// @param yrhole    r-dependent embedded DOS [251][lx+1] (flat)
/// @param yrhoce    r-dependent central DOS [251] (accumulated)
/// @param yrhocp    r-dependent previous DOS [251] (accumulated)
/// @param rhoval    Valence density [251] (output, accumulated)
/// @param xnmues    Occupation numbers [lx+1] (accumulated)
/// @param xnatph    Number of atoms for this potential
/// @param xntot     Total occupation (accumulated, output)
/// @param iflr      Flag: 1 = add contribution from point to real axis
/// @param iflrp     Flag: 1 = add contribution from previous point to real axis
/// @param fl        Running integral left (accumulated)
/// @param fr        Running integral right (accumulated)
/// @param iunf      Unfreeze f-electrons flag
void ff2g(const std::complex<float>* gtr, int iph, int ie, int ilast,
          FeffComplex* xrhoce, const FeffComplex* xrhole,
          FeffComplex* xrhocp,
          FeffComplex ee, FeffComplex ep,
          const FeffComplex* yrhole, FeffComplex* yrhoce, FeffComplex* yrhocp,
          double* rhoval, double* xnmues,
          double xnatph, double& xntot, int iflr, int iflrp,
          FeffComplex& fl, FeffComplex& fr, int iunf);

} // namespace feff::pot
