#pragma once
// L-resolved DOS and phase shifts via Dirac equation.
// Converted from src/POT/rholie.f
//
// Computes the l-resolved density of states and phase shifts by solving
// the Dirac equation twice per l-channel: once for the regular solution
// and once for the irregular solution. Calls dfovrg() from FOVRG module.

#include <feff/dimensions.hpp>
#include <feff/types.hpp>

namespace feff::pot {

/// Compute l-resolved DOS and phase shifts at complex energy em.
///
/// @param ri05     Radial grid with 0.05 step [251] (output: filled)
/// @param nr05     Output: number of valid points in ri05
/// @param dx       Grid spacing for Loucks grid
/// @param x0       Grid origin parameter (8.8)
/// @param ri       Loucks radial grid [nrptx]
/// @param em       Complex energy (Hartrees)
/// @param ixc      Exchange-correlation model
/// @param rmt      Muffin-tin radius
/// @param rnrm     Norman radius
/// @param vtot     Total potential on Loucks grid [nrptx] (complex, will be cast)
/// @param vvalgs   Valence potential on Loucks grid [nrptx] (complex, will be cast)
/// @param xnval    Valence occupation numbers [30]
/// @param dgcn     Large Dirac components on Loucks grid [nrptx][30] (flat)
/// @param dpcn     Small Dirac components on Loucks grid [nrptx][30] (flat)
/// @param eref     Reference energy (complex)
/// @param adgc     Development coefficients [10][30] (flat)
/// @param adpc     Development coefficients [10][30] (flat)
/// @param xrhole   Output: integrated l-DOS for embedded atom [lx+1]
/// @param xrhoce   Output: integrated l-DOS for central atom [lx+1]
/// @param yrhole   Output: r-dependent l-DOS [251][lx+1] (flat)
/// @param yrhoce   Output: r-dependent central DOS [251]
/// @param ph       Output: phase shifts [lx+1]
/// @param iz       Atomic number
/// @param xion     Ionicity
/// @param iunf     Unfreeze f-electrons flag
/// @param ihole    Core-hole orbital index
/// @param lmaxsc   Maximum angular momentum for scattering
void rholie(double* ri05, int& nr05, double dx, double x0,
            double* ri, FeffComplex em,
            int ixc, double rmt, double rnrm,
            double* vtot, double* vvalgs,
            const double* xnval, double* dgcn, double* dpcn,
            FeffComplex eref,
            const double* adgc, const double* adpc,
            FeffComplex* xrhole, FeffComplex* xrhoce,
            FeffComplex* yrhole, FeffComplex* yrhoce,
            FeffComplex* ph,
            int iz, double xion, int iunf, int ihole, int lmaxsc);

} // namespace feff::pot
