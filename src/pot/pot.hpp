#pragma once
// Main POT computational kernel.
// Converted from: src/POT/pot.f (~680 lines)
//
// Computes potentials for an input atomic cluster, returning data needed
// to compute phase shifts. Performs:
//   - Free atom calculations (scfdat)
//   - Overlap potentials/densities (ovrlp)
//   - Muffin-tin radii and interstitial parameters (istprm, afolp)
//   - Core-valence separation (corval)
//   - Self-consistent field loop (scmt/scmtmp)
//   - Optional diagnostic output (wpot)

#include <feff/dimensions.hpp>

namespace feff::pot {

/// Main POT computational kernel.
///
/// Array layout (Fortran column-major, flat pointers):
///   rat[3*natx]                    - atomic coordinates (column-major: x,y,z per atom)
///   iphat[natx]                    - potential index per atom (0-based pot indices)
///   iatph[nphx+1]                  - representative atom for each potential
///   xnatph[nphx+1]                 - stoichiometry per potential type
///   iphovr[novrx*(nphx+1)]         - overlap potential indices
///   nnovr[novrx*(nphx+1)]          - overlap atom counts
///   rovr[novrx*(nphx+1)]           - overlap radii
///   novr[nphx+1]                   - number of overlap shells per potential
///   folp0[nphx+1], xion[nphx+1]   - overlap factors, ionicities
///   iz[nphx+1]                     - atomic numbers
///   lmaxsc[nphx+1]                 - max l for SCF per potential
///   imt[nphx+1], rmt[nphx+1]      - MT mesh index and radius (output)
///   inrm[nphx+1], rnrm[nphx+1]    - Norman mesh index and radius (output)
///   folpx[nphx+1]                  - max overlap factor (output)
///   dgc0[251], dpc0[251]           - initial orbital spinors (output)
///   dgc[251*30*(nphx+2)]           - Dirac upper spinor components (output)
///   dpc[251*30*(nphx+2)]           - Dirac lower spinor components
///   adgc[10*30*(nphx+2)]           - development coefficients
///   adpc[10*30*(nphx+2)]           - development coefficients
///   edens[251*(nphx+1)]            - total electron density (output)
///   vclap[251*(nphx+1)]            - overlap Coulomb potential (output)
///   vtot[251*(nphx+1)]             - total potential (output)
///   edenvl[251*(nphx+1)]           - valence density (output)
///   vvalgs[251*(nphx+1)]           - valence GS potential (output)
///   dmag[251*(nphx+2)]             - spin magnetization density
///   xnval[30*(nphx+2)]             - valence electron counts
///   eorb[30*(nphx+2)]              - orbital energies (output)
///   kappa[30*(nphx+2)]             - kappa quantum numbers (output)
///   iorb[8*(nphx+2)]               - orbital index (-4:3 per iph, stride 8)
///   qnrm[nphx+1]                   - charge within Norman sphere (output)
///   xnmues[(lx+1)*(nphx+1)]        - SCF occupation numbers (output)
///   title[nheadx]                   - title lines (80 chars each)
void pot(bool verbose, double rgrd, int nohole,
         int inters, double totvol, double ecv0, int nscmt, int nmix,
         int ntitle, char title[][80],
         int nat, int nph, int ihole, int iafolp, int ixc,
         int* iphat, double* rat, int* iatph, double* xnatph,
         int* novr, int* iphovr, int* nnovr, double* rovr,
         double* folp0, double* xion, int iunf, int* iz, int ipr1,
         int ispec, int jumprm, int* lmaxsc, int icoul,
         double ca1, float rfms1, int lfms1,
         // output scalars
         double& rnrmav, double& xmu, double& vint, double& rhoint,
         double& emu, double& s02, double& erelax, double& wp,
         double& rs, double& xf, double& qtotel,
         // output arrays
         int* imt, double* rmt, int* inrm, double* rnrm, double* folpx,
         double* dgc0, double* dpc0,
         double* dgc, double* dpc, double* adgc, double* adpc,
         double* edens, double* vclap, double* vtot,
         double* edenvl, double* vvalgs, double* dmag, double* xnval,
         double* eorb, int* kappa, int* iorb,
         double* qnrm, double* xnmues, int& nhtmp);

} // namespace feff::pot
