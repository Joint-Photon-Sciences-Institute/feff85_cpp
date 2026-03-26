#pragma once
// Serial self-consistency loop for muffin-tin potentials.
// Converted from: src/POT/scmt.f
//
// Finds new Fermi level (xmu), electron counts (qnrm), and
// new valence densities (rhoval) via complex energy contour integration.

#include <feff/dimensions.hpp>

namespace feff::pot {

/// Serial SCF iteration in the complex energy plane.
///
/// Array layout (Fortran-order, flat pointers):
///   edens[251*(nphx+1)]       - electron density, stride 251 per iph
///   edenvl[251*(nphx+1)]      - valence density
///   vtot[251*(nphx+1)]        - total potential
///   vvalgs[251*(nphx+1)]      - valence GS potential
///   vclap[251*(nphx+1)]       - overlap Coulomb potential
///   rhoval[251*(nphx+2)]      - valence density (output), stride 251 per iph
///   dgc[251*30*(nphx+2)]      - Dirac upper spinor, stride 251*30 per iph
///   dpc[251*30*(nphx+2)]      - Dirac lower spinor
///   adgc[10*30*(nphx+2)]      - development coeff upper
///   adpc[10*30*(nphx+2)]      - development coeff lower
///   xnval[30*(nphx+2)]        - valence electron counts
///   xnvmu[(lx+1)*(nphx+2)]    - occupation numbers from getorb
///   xnmues[(lx+1)*(nphx+1)]   - occupation numbers from SCF (output)
///   rat[3*natx]                - atomic coordinates
///   lmaxsc[nphx+1]            - max l for SCF per potential
///
/// Scalars (0-based potential arrays where noted):
///   rmt[nphx+1], rnrm[nphx+1], qnrm[nphx+1] - MT radii, Norman radii, charges
///   xnatph[nphx+1], xion[nphx+1], iz[nphx+1] - stoichiometry, ionicity, atomic number
///   iatph[nphx+1]                              - representative atom for each potential
///   iphat[natx]                                 - potential index per atom
void scmt(bool verbose, int iscmt, double ecv, int nph, int nat,
          double* vclap, double* edens, double* edenvl,
          double* vtot, double* vvalgs,
          double* rmt, double* rnrm, double* qnrm,
          int ixc, double rhoint, double vint, double& xmu, int jumprm,
          double xnferm, double* xnvmu, double* xnval,
          double x0, double* ri, double dx,
          double* xnatph, double* xion, int iunf, int* iz,
          double* adgc, double* adpc, double* dgc, double* dpc,
          int ihole,
          double* rat, int* iatph, int* iphat,
          int* lmaxsc, double* rhoval, double* xnmues,
          bool& ok,
          double rgrd, int nohole, int nscmt, int icoul,
          double ca1, float rfms1, int lfms1);

} // namespace feff::pot
