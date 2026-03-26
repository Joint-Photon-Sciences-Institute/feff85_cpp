#pragma once
// Read input files for POT module.
// Converted from src/POT/reapot.f
//
// Reads geometry (geom.json) and potential parameters (pot.json)
// using JSON readers, then converts to code units (Bohr and Hartrees).

#include <feff/dimensions.hpp>
#include <string>

namespace feff::pot {

/// Read all POT module input from JSON files and convert to atomic units.
///
/// @param mpot     Potential generation flag (output)
/// @param rgrd     Grid spacing (output)
/// @param ntitle   Number of title lines (output)
/// @param title    Title lines [nheadx] (output)
/// @param ipr1     Print flag (output)
/// @param ispec    Spectroscopy type (output)
/// @param nohole   Core-hole flag (output)
/// @param ihole    Core-hole orbital index (output)
/// @param gamach   Core-hole lifetime broadening [eV -> Hartrees] (output)
/// @param nph      Number of unique potentials (output)
/// @param iz       Atomic numbers [nphx+1] (output)
/// @param lmaxsc   Max angular momentum [nphx+1] (output)
/// @param xnatph   Atoms per potential type [nphx+1] (output)
/// @param xion     Ionicity [nphx+1] (output)
/// @param iunf     Unfreeze f-electrons flag (output)
/// @param ixc      Exchange-correlation model (output)
/// @param jumprm   Jump flag (output)
/// @param iafolp   Auto-FOLP flag (output)
/// @param folp     Overlap fractions [nphx+1] (output)
/// @param inters   Interstitial model (output)
/// @param totvol   Total volume [Ang^3 -> Bohr^3] (output)
/// @param rfms1    FMS cluster radius [Ang -> Bohr] (output)
/// @param lfms1    FMS l-max (output)
/// @param nscmt    Number of SCF iterations (output)
/// @param ca1      Convergence accelerator (output)
/// @param nmix     Mixing parameter (output)
/// @param ecv      Core-valence separation [eV -> Hartrees] (output)
/// @param icoul    Coulomb normalization mode (output)
/// @param novr     Overlap neighbor counts [nphx+1] (output)
/// @param iphovr   Overlap neighbor pot types [novrx][nphx+1] (flat, output)
/// @param nnovr    Overlap neighbor multiplicities [novrx][nphx+1] (flat, output)
/// @param rovr     Overlap neighbor distances [novrx][nphx+1] (flat, output)
/// @param nat      Number of atoms (output)
/// @param rat      Atomic coordinates [3][natx] (column-major, Ang -> Bohr, output)
/// @param iphat    Potential type per atom [natx] (output)
/// @param iatph    Representative atom per potential [nphx+1] (output)
void reapot(int& mpot, double& rgrd, int& ntitle, std::string* title,
            int& ipr1, int& ispec, int& nohole, int& ihole, double& gamach,
            int& nph, int* iz, int* lmaxsc, double* xnatph,
            double* xion, int& iunf, int& ixc, int& jumprm, int& iafolp,
            double* folp, int& inters, double& totvol,
            float& rfms1, int& lfms1, int& nscmt, double& ca1, int& nmix,
            double& ecv, int& icoul,
            int* novr, int* iphovr, int* nnovr, double* rovr,
            int& nat, double* rat, int* iphat, int* iatph);

} // namespace feff::pot
