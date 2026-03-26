#pragma once
// Find MT radii and interstitial parameters.
// Converted from src/POT/istprm.f
//
// Determines muffin-tin radii via Norman prescription, computes
// interstitial density and potential, and finds the Fermi level.
// Calls fermi(), vbh(), edp(), movrlp(), ovp2mt().

#include <feff/dimensions.hpp>
#include <feff/types.hpp>

namespace feff::pot {

/// Compute interstitial parameters, MT radii, Fermi level, etc.
///
/// @param nph      Number of unique potentials (0-based max index)
/// @param nat      Number of atoms in cluster
/// @param iphat    Potential type per atom [natx]
/// @param rat      Atomic coordinates [3][natx] (column-major)
/// @param iatph    Representative atom for each potential [nphx+1]
/// @param xnatph   Number of atoms per potential type [nphx+1]
/// @param novr     Number of explicit overlap neighbors [nphx+1]
/// @param iphovr   Pot type of overlap neighbors [novrx][nphx+1] (flat)
/// @param nnovr    Count of overlap neighbors [novrx][nphx+1] (flat)
/// @param rovr     Distance to overlap neighbors [novrx][nphx+1] (flat)
/// @param folp     Overlap fractions [nphx+1]
/// @param folpx    Maximum overlap fractions [nphx+1] (output)
/// @param iafolp   Automatic FOLP flag
/// @param edens    Overlapped electron density [251][nphx+1] (flat)
/// @param edenvl   Overlapped valence density [251][nphx+1] (flat)
/// @param idmag    Spin-dependent flag
/// @param dmag     Spin density ratio [251][nphx+2] (flat)
/// @param vclap    Overlapped Coulomb potential [251][nphx+1] (flat)
/// @param vtot     Total potential (output) [251][nphx+1] (flat)
/// @param vvalgs   Valence gs potential (output) [251][nphx+1] (flat)
/// @param imt      MT grid indices [nphx+1] (output)
/// @param inrm     Norman grid indices [nphx+1] (output)
/// @param rmt      MT radii [nphx+1] (input/output)
/// @param rnrm     Norman radii [nphx+1] (input/output)
/// @param ixc      Exchange-correlation model
/// @param rhoint   Interstitial density (output)
/// @param vint     Interstitial potential (output)
/// @param rs       Density parameter (output)
/// @param xf       Fermi momentum (output)
/// @param xmu      Fermi level (input, previous value)
/// @param xmunew   New Fermi level (output)
/// @param rnrmav   Average Norman radius (output)
/// @param qtotel   Total cluster charge (output)
/// @param inters   Interstitial model flag (input/output)
/// @param totvol   Total volume (input)
void istprm(int nph, int nat, const int* iphat, const double* rat,
            const int* iatph, const double* xnatph,
            const int* novr, const int* iphovr, const int* nnovr, const double* rovr,
            const double* folp, double* folpx, int iafolp,
            double* edens, double* edenvl, int idmag,
            double* dmag, double* vclap, double* vtot, double* vvalgs,
            int* imt, int* inrm, double* rmt, double* rnrm,
            int ixc, double& rhoint, double& vint, double& rs, double& xf,
            double xmu, double& xmunew,
            double& rnrmav, double& qtotel, int& inters, double totvol);

/// Calculate lens volume for overlapping spheres (local helper).
/// vol_i = (pi/3) * h^2 * (3*r1 - h)  where h = r1 - xl
/// and xl = (r1^2 - r2^2 + r^2) / (2*r)
double calcvl(double r1, double r2, double r);

} // namespace feff::pot
