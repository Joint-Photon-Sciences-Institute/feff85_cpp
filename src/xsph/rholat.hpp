#pragma once
// Spectroscopic density of states (Renormalized atom method).
// Converted from src/XSPH/rholat.f

#include <feff/types.hpp>
#include <feff/dimensions.hpp>

namespace feff::xsph {

/// Calculate spectroscopic DOS using renormalized atom method.
///
/// @param icount  1=renormalized atom counts, 2=Mulliken counts
/// @param dx      Grid step
/// @param x0      Grid parameter
/// @param ri      Radial grid [nrptx]
/// @param em      Complex energy point
/// @param ixc     Exchange-correlation model
/// @param rmt     Muffin-tin radius
/// @param rnrm    Norman radius
/// @param vtot    Total potential [nrptx]
/// @param vvalgs  Valence potential [nrptx]
/// @param xnval   Valence occupations [30]
/// @param iorb    Orbital indices [-4:3] (8 elements, offset by 4)
/// @param dgcn    Large Dirac components [nrptx][30]
/// @param dpcn    Small Dirac components [nrptx][30]
/// @param eref    Energy reference (complex, in/out)
/// @param adgc    Development coefficients [10][30]
/// @param adpc    Development coefficients [10][30]
/// @param xrhole  Output: density integral [-4:3][-4:3] (8x8 complex)
/// @param xrhoce  Output: central atom density [-4:3][-4:3] (8x8 complex)
/// @param ph      Output: phase shifts [lx+1]
/// @param iz      Atomic number
/// @param xion    Ionicity
/// @param iunf    Unfolding flag
/// @param ihole   Core-hole index
/// @param lmaxsc  Max l for scattering
void rholat(int icount, double dx, double x0, const double ri[],
            FeffComplex em, int ixc, double rmt, double rnrm,
            const double vtot[], const double vvalgs[],
            const double xnval[], const int iorb[],
            const double dgcn[][30], const double dpcn[][30],
            FeffComplex eref,
            const double adgc[][30], const double adpc[][30],
            FeffComplex* xrhole, FeffComplex* xrhoce,
            FeffComplex ph[], int iz, double xion, int iunf,
            int ihole, int lmaxsc);

} // namespace feff::xsph
