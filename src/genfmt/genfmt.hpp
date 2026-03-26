#pragma once

// Main genfmt F-matrix generator.
// Converted from GENFMT/genfmt.f (~483 lines)
//
// Computes EXAFS F-matrix scattering amplitudes for all paths,
// writes feff.pad and list.dat output files.

#include "genfmt_data.hpp"
#include "regenf.hpp"

namespace feff::genfmt {

/// Main genfmt routine: compute F-matrix for all paths.
///
/// ipr5:   print level for feff.dat output
/// critcw: curved-wave chi amplitude ratio threshold (%)
/// iorder: order of approximation in f-matrix expansion
/// wnstar: flag to write nstar.dat
/// ipol:   polarization type
/// ispin:  spin flag
/// le2:    multipole moment selector
/// angks:  angle between k-vector and spin-vector
/// elpty:  ellipticity
/// evec:   polarization vector
/// xivec:  direction of travel
/// ptz:    polarization tensor
void genfmt(int ipr5, double critcw, int iorder, bool wnstar,
            int ipol, int ispin, int le2, double angks, double elpty,
            const double evec[3], const double xivec[3],
            const FeffComplex ptz[3][3]);

} // namespace feff::genfmt
