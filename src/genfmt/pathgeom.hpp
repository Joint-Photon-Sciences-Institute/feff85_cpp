#pragma once

// Path geometry computation (Euler angles and leg lengths).
// Converted from GENFMT/pathgeom.f
//
// Using the coordinates of the path's constituent atoms,
// computes ri (leg lengths), beta (Euler beta angles), and eta (combined Euler angles).

#include <feff/dimensions.hpp>

namespace feff::genfmt {

/// Compute path geometry from atom positions.
///
/// Input:
///   nleg: number of legs in the path
///   ipol: polarization flag (> 0 for polarization case)
///   rat: atom positions rat[xyz][atom_idx], atom_idx in [0, legtot+1]
///   ipot: potential indices ipot[0..legtot]
///
/// Output:
///   nsc: number of scatterers (nleg - 1)
///   ri: leg lengths ri[0..nleg-1]
///   beta: Euler beta angles beta[0..nleg] (or nleg+1 for polarization)
///   eta: combined Euler angles eta[0..legtot+1]
void pathgeom(int nleg, int& nsc, int ipol,
              double rat[][legtot + 2], int ipot[],
              double ri[], double beta[], double eta[]);

} // namespace feff::genfmt
