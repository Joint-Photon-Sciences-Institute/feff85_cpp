#pragma once
// Core-valence separation.
// Converted from src/POT/corval.f
//
// Finds the core-valence separation energy (ecv) for the cluster.
// Searches for suspicious peaks in the l-resolved DOS between -20 and -70 eV,
// sorts them, and determines which states are core vs valence.
// Written by ALA 10/1998.

#include <feff/dimensions.hpp>
#include <feff/types.hpp>

namespace feff::pot {

/// Find core-valence separation for the cluster of atoms.
///
/// @param verbse   Verbose output flag
/// @param ecv      Core-valence separation energy (input/output) [Hartrees]
/// @param xnvmu    Valence electron counts [(lx+1)][nphx+2] (flat, output)
/// @param eorb     Orbital energies [30][nphx+1] (flat)
/// @param norb     Number of orbitals per potential [nphx+1]
/// @param xnval    Valence occupation numbers [30][nphx+1] (flat, modified)
/// @param kappa    Kappa quantum numbers [30][nphx+1] (flat)
/// @param rgrd     Grid spacing for Loucks grid
/// @param nohole   Core-hole flag
/// @param nph      Number of unique potentials
/// @param edens    Electron density [251][nphx+1] (flat)
/// @param edenvl   Valence density [251][nphx+1] (flat, modified)
/// @param vtot     Total potential [251][nphx+1] (flat)
/// @param vvalgs   Valence potential [251][nphx+1] (flat)
/// @param rmt      MT radii [nphx+1]
/// @param rnrm     Norman radii [nphx+1]
/// @param ixc      Exchange-correlation model
/// @param rhoint   Interstitial density
/// @param vint     Interstitial potential
/// @param jumprm   Jump flag
/// @param x0       Grid origin parameter (8.8)
/// @param ri       Radial grid work array [nrptx]
/// @param dx       Grid spacing
/// @param xion     Ionicity per potential [nphx+1]
/// @param iunf     Unfreeze f-electrons flag
/// @param iz       Atomic numbers [nphx+1]
/// @param adgc     Development coefficients for large component [10][30][nphx+2] (flat)
/// @param adpc     Development coefficients for small component [10][30][nphx+2] (flat)
/// @param dgc      Large Dirac component [251][30][nphx+2] (flat)
/// @param dpc      Small Dirac component [251][30][nphx+2] (flat)
/// @param ihole    Core-hole orbital index
/// @param lmaxsc   Maximum angular momentum per potential [nphx+1]
void corval(bool verbse,
            double& ecv, double* xnvmu, const double* eorb, const int* norb,
            double* xnval, const int* kappa, double rgrd,
            int nohole, int nph, double* edens, double* edenvl,
            const double* vtot, const double* vvalgs,
            const double* rmt, const double* rnrm, int ixc,
            double rhoint, double vint, int jumprm,
            double x0, double* ri, double dx, const double* xion,
            int iunf, const int* iz,
            const double* adgc, const double* adpc,
            double* dgc, double* dpc,
            int ihole, const int* lmaxsc);

} // namespace feff::pot
