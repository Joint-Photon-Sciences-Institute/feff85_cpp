#pragma once
// Grid interpolation utilities for Dirac spinors and potentials.
// Converted from: src/COMMON/fixdsp.f, src/COMMON/fixdsx.f, src/COMMON/fixvar.f

#include <feff/dimensions.hpp>

namespace feff::common {

/// Interpolate Dirac spinor components (dgc, dpc) from ATOM grid to xsect grid.
/// Converted from fixdsp.f.
///   dxorg, dxnew: grid spacings (original and new)
///   dgc0[251], dpc0[251]: upper/lower spinor components on original grid
///   dgcx[nrptx], dpcx[nrptx]: interpolated components on new grid (output)
///   jnew: number of valid points in new grid (output)
void fixdsp(double dxorg, double dxnew,
            const double dgc0[251], const double dpc0[251],
            double dgcx[], double dpcx[], int& jnew);

/// Interpolate Dirac spinor components for multiple orbitals for one potential.
/// Converted from fixdsx.f.
///   iph: potential index
///   dxorg, dxnew: grid spacings
///   dgc[251][30][nphx+2], dpc[251][30][nphx+2]: all spinor components
///   dgcn[nrptx][30], dpcn[nrptx][30]: interpolated components (output)
void fixdsx(int iph, double dxorg, double dxnew,
            const double* dgc, const double* dpc,
            double* dgcn, double* dpcn);

/// Interpolate potentials and electron densities to new grid for PHASE code.
/// Converted from fixvar.f.
///   rmt: muffin-tin radius
///   edens[251], vtot[251], dmag[251]: density, potential, magnetization on original grid
///   vint, rhoint: interstitial potential and density
///   dxorg, dxnew: grid spacings
///   jumprm: potential jump flag (1=calculate vjump, >0=apply vjump)
///   vjump: potential discontinuity at MT boundary (input/output)
///   ri[nrptx]: radial grid (output)
///   vtotph[nrptx], rhoph[nrptx], dmagx[nrptx]: interpolated arrays (output)
void fixvar(double rmt, const double edens[251], const double vtot[251],
            const double dmag[251], double vint, double rhoint,
            double dxorg, double dxnew, int jumprm,
            double& vjump, double ri[], double vtotph[],
            double rhoph[], double dmagx[]);

} // namespace feff::common
