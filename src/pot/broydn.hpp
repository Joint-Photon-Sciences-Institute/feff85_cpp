#pragma once
// Broyden density mixing algorithm.
// Converted from src/POT/broydn.f
//
// Calculates new density using Broyden algorithm
// (J.Phys.A, 17, L317, 1984).
// Also handles the charge inside each Norman sphere properly.

#include <feff/dimensions.hpp>

namespace feff::pot {

/// Broyden density mixing for SCF convergence.
///
/// @param iscmt    SCF iteration number (1-based)
/// @param ca       Convergence accelerator factor
/// @param nph      Number of unique potentials
/// @param xnvmu    Valence electron counts [(lx+1)][nphx+2] (flat)
/// @param ilast    Last grid index per potential [nphx+1]
/// @param xnatph   Atoms per potential type [nphx+1]
/// @param rnrm     Norman radii [nphx+1]
/// @param qnrm     Charge inside Norman sphere [nphx+1] (modified)
/// @param edenvl   Old valence density [251][nphx+1] (flat)
/// @param rhoval   New valence density from integration [251][nphx+2] (flat, modified)
///                 At input: density*4*pi*r^2; at output: mixed density
/// @param dq       Charge transfer per potential [nphx+1] (output)
void broydn(int iscmt, double ca, int nph, const double* xnvmu,
            const int* ilast, const double* xnatph, const double* rnrm,
            double* qnrm, double* edenvl, double* rhoval, double* dq);

} // namespace feff::pot
