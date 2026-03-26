#pragma once
// Cross-section calculation (absorption and scattering).
// Converted from src/XSPH/xsect.f (~761 lines)
//
// Calculates atomic absorption cross-section and reduced dipole/multipole
// matrix elements for the absorbing atom.

#include <feff/types.hpp>
#include <feff/dimensions.hpp>

namespace feff::xsph {

/// Calculate cross-section and reduced matrix elements.
///
/// @param ipr2    Print level
/// @param dx      Loucks grid step
/// @param x0      Grid parameter
/// @param ri      Radial grid [nrptx]
/// @param ne      Total energy points
/// @param ne1     Horizontal grid points
/// @param ik0     Index where k=0
/// @param em      Complex energy grid [nex]
/// @param edge    Chemical potential (Hartrees)
/// @param ihole   Core-hole index
/// @param emu     Edge energy (Hartrees)
/// @param corr    Correction factor
/// @param dgc0    Large Dirac component, initial orbital [nrptx]
/// @param dpc0    Small Dirac component, initial orbital [nrptx]
/// @param jnew    Number of valid points in interpolated initial orbital
/// @param ixc     Exchange-correlation model
/// @param lreal   Real phase shift flag
/// @param rmt     Muffin-tin radius
/// @param rnrm    Norman radius
/// @param xmu     Fermi level (Hartrees)
/// @param iPl     Many-pole control
/// @param vtot    Total potential [nrptx]
/// @param vvalgs  Valence potential [nrptx]
/// @param edens   Electron density [nrptx]
/// @param dmag    Magnetization density [nrptx]
/// @param edenvl  Valence density [nrptx]
/// @param dgcn    Large Dirac components [nrptx][30]
/// @param dpcn    Small Dirac components [nrptx][30]
/// @param adgc    Development coefficients [10][30]
/// @param adpc    Development coefficients [10][30]
/// @param xsec    Output: atomic cross-section [nex]
/// @param xsnorm  Output: normalization cross-section [nex]
/// @param rkk     Output: reduced matrix elements [nex][8]
/// @param iz      Atomic number
/// @param xion    Ionicity
/// @param iunf    Unfolding flag
/// @param xnval   Valence occupations [30]
/// @param izstd   TDLDA flag
/// @param iorb    Orbital indices [-4:3] (8 elements)
/// @param l2lp    l-->l+1 / l-->l-1 transition selector
/// @param ipol    Polarization type
/// @param ispin   Spin flag
/// @param le2     Multipole selector (0=E1, 1=M1, 2=E2)
/// @param angks   k-spin angle
/// @param ptz     Polarization tensor [3][3]
void xsect(int ipr2, double dx, double x0, const double ri[],
           int ne, int ne1, int ik0, const FeffComplex em[], double edge,
           int ihole, double emu, double corr,
           const double dgc0[], const double dpc0[], int jnew,
           int ixc, int lreal, double rmt, double rnrm, double xmu,
           int iPl,
           const double vtot[], const double vvalgs[],
           const double edens[], const double dmag[], const double edenvl[],
           const double dgcn[][30], const double dpcn[][30],
           const double adgc[][30], const double adpc[][30],
           FeffComplex xsec[], double xsnorm[], FeffComplex rkk[][8],
           int iz, double xion, int iunf, const double xnval[],
           int izstd, const int iorb[], int l2lp,
           int ipol, int ispin, int le2, double angks,
           const FeffComplex ptz[3][3]);

} // namespace feff::xsph
