#pragma once
// Phase shift calculation per unique potential.
// Converted from src/XSPH/phase.f
//
// Calculates complex scattering phase shifts for a given potential
// by solving the Dirac equation and matching to free-particle solutions
// at the muffin-tin boundary.

#include <feff/types.hpp>
#include <feff/dimensions.hpp>

namespace feff::xsph {

/// Calculate phase shifts for one unique potential.
///
/// @param iph     Unique potential index (for messages only)
/// @param dx      Loucks grid step
/// @param x0      Grid parameter (ri = exp((i-1)*dx - x0))
/// @param ri      Radial grid [nrptx]
/// @param ne      Total energy points
/// @param ne1     Horizontal grid points
/// @param ne3     Auxiliary horizontal points
/// @param em      Complex energy grid [nex]
/// @param ixc     Exchange-correlation model
/// @param nsp     Number of spin channels
/// @param lreal   Real self-energy flag
/// @param rmt     Muffin-tin radius
/// @param rnrm    Norman radius
/// @param xmu     Fermi level (Hartrees)
/// @param iPl     Many-pole self-energy control
/// @param vtot    Total potential [nrptx]
/// @param vvalgs  Valence potential [nrptx]
/// @param edens   Electron density [nrptx]
/// @param dmag    Magnetization density [nrptx]
/// @param edenvl  Valence density [nrptx]
/// @param dgcn    Large Dirac components [nrptx][30]
/// @param dpcn    Small Dirac components [nrptx][30]
/// @param adgc    Development coefficients [10][30]
/// @param adpc    Development coefficients [10][30]
/// @param eref    Output: energy reference [nex]
/// @param ph      Output: phase shifts [nex][-ltot:ltot]
/// @param lmax    Output: max angular momentum used
/// @param iz      Atomic number
/// @param ihole   Core-hole index (nonzero for absorber)
/// @param xion    Ionicity
/// @param iunf    Unfolding flag
/// @param xnval   Valence occupations [30]
/// @param ispin   Spin flag
void phase(int iph, double dx, double x0, const double ri[],
           int ne, int ne1, int ne3, const FeffComplex em[],
           int ixc, int nsp, int lreal, double rmt, double rnrm,
           double xmu, int iPl,
           const double vtot[], const double vvalgs[],
           const double edens[], const double dmag[], const double edenvl[],
           const double dgcn[][30], const double dpcn[][30],
           const double adgc[][30], const double adpc[][30],
           FeffComplex eref[], FeffComplex* ph, int& lmax,
           int iz, int ihole, double xion, int iunf,
           const double xnval[], int ispin);

} // namespace feff::xsph
