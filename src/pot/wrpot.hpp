#pragma once
// Write pot.pad output file in PAD (Packed ASCII Data) format.
// Converted from: src/POT/wrpot.f
//
// Writes all potential data needed by downstream modules (XSPH, etc.)
// to a portable ASCII file.

#include <feff/dimensions.hpp>
#include <string>

namespace feff::pot {

/// Write pot.pad file containing all potential/density data.
///
/// Array layout (Fortran column-major, flat pointers):
///   imt[nphx+1], rmt[nphx+1]       - MT mesh index and radius
///   inrm[nphx+1], rnrm[nphx+1]     - Norman mesh index and radius
///   folp[nphx+1], folpx[nphx+1]    - overlap factors
///   xnatph[nphx+1]                   - stoichiometry
///   xion[nphx+1], iz[nphx+1]        - ionicities, atomic numbers
///   dgc0[251], dpc0[251]             - initial orbital spinors
///   dgc[251*30*(nphx+1)]             - Dirac upper spinor (0:nphx only)
///   dpc[251*30*(nphx+1)]             - Dirac lower spinor
///   adgc[10*30*(nphx+1)]             - development coefficients
///   adpc[10*30*(nphx+1)]             - development coefficients
///   edens[251*(nphx+1)]              - total electron density
///   vclap[251*(nphx+1)]              - Coulomb potential
///   vtot[251*(nphx+1)]               - total potential
///   edenvl[251*(nphx+1)]             - valence density
///   vvalgs[251*(nphx+1)]             - valence GS potential
///   dmag[251*(nphx+1)]               - spin magnetization
///   xnval[30*(nphx+1)]               - valence electron counts
///   eorb[30]                          - orbital energies (absorber only)
///   kappa[30]                         - kappa quantum numbers (absorber only)
///   iorb[8*(nphx+1)]                 - orbital indices (-4:3 per iph)
///   qnrm[nphx+1]                     - charges within Norman sphere
///   xnmues[(lx+1)*(nphx+1)]          - SCF occupation numbers
///   title[nheadx]                     - title lines (80 chars each)
void wrpot(int nph, int ntitle, const char title[][80],
           double rnrmav, double xmu, double vint, double rhoint,
           double emu, double s02, double erelax, double wp,
           double ecv, double rs, double xf, double qtotel,
           const int* imt, const double* rmt, const int* inrm, const double* rnrm,
           const double* folp, const double* folpx, const double* xnatph,
           const double* dgc0, const double* dpc0,
           const double* dgc, const double* dpc,
           const double* adgc, const double* adpc,
           const double* edens, const double* vclap, const double* vtot,
           const double* edenvl, const double* vvalgs, const double* dmag,
           const double* xnval,
           const double* eorb, const int* kappa, const int* iorb,
           const double* qnrm, const double* xnmues,
           int nohole, int ihole,
           int inters, double totvol, int iafolp,
           const double* xion, int iunf, const int* iz, int jumprm);

} // namespace feff::pot
