#pragma once
// LDA exchange-correlation potential with core-valence separation.
// Converted from: src/ATOM/vlda.f
//
// Note: Calls vbh() and edp() from EXCH module. Currently uses stub
// implementations that will be replaced when EXCH is converted (Phase 4).

#include "atom_types.hpp"

namespace feff::atom {

/// Calculate LDA exchange-correlation potential and add to direct potential.
/// Computes total and valence electron densities, then applies xc model.
///   xnval[30]: valence occupation numbers
///   srho[251]: total electron density (output)
///   srhovl[251]: valence electron density (output)
///   vtrho[251]: accumulated V_xc * rho for total energy (input/output)
///   ilast: if > 0, accumulate energy contribution
///   idfock: exchange model (1=DF, 2=LDA, 5=exch5, 6=exch6)
/// Replaces Fortran vlda().
void vlda(const double xnval[30], double srho[251], double srhovl[251],
          double vtrho[251], int ilast, int idfock, AtomState& state);

// =========================================================================
// EXCH module stubs — will be replaced in Phase 4
// =========================================================================

/// Vosko-Wilk-Nusair exchange-correlation potential.
/// rs: Wigner-Seitz radius, xm: spin polarization (usually 1.0)
/// vxc: exchange-correlation potential in Hartrees (output)
void vbh_stub(double rs, double xm, double& vxc);

/// Dirac-Hara exchange potential correction.
/// rs: Wigner-Seitz radius, xf: Fermi momentum
/// vdh: Dirac-Hara exchange potential in Hartrees (output)
void edp_stub(double rs, double xf, double& vdh);

} // namespace feff::atom
