#pragma once
// Make path parameters xp, yp, zp for each atom in standard frame.
// Converted from: src/PATH/mpprmp.f
//
// Computes coordinates of path atoms in a standard reference frame
// determined by the polarization case (icase).

#include "path_data.hpp"

namespace feff::path {

/// Compute x,y,z coordinates for a path in the standard frame of reference.
/// npat: number of path atoms.
/// ipat[0..npat]: atom indices (ipat[npat] = 0 for central atom).
/// xp, yp, zp[0..npat-1]: output coordinates.
/// ipol, ispin: polarization and spin flags.
/// evec[3], xivec[3]: polarization and incidence vectors.
/// ica: override for icase if positive (for EELS).
void mpprmp(int npat, const int ipat[], float xp[], float yp[], float zp[],
            int ipol, int ispin, const double evec[3], const double xivec[3],
            int ica, const AtomData& atoms);

} // namespace feff::path
