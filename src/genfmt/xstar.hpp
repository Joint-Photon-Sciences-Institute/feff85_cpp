#pragma once

// Dichroism / n* calculation.
// Converted from GENFMT/xstar.f
//
// Calculates nstar = deg * cos(eps,r1) * cos(eps,rN) using plane wave
// approximation for the central atom.
// Written by Alexei Ankudinov, 08.13.96

#include <feff/dimensions.hpp>

namespace feff::genfmt {

/// Calculate the n* factor for a path.
///
/// eps1: polarization vector
/// eps2: ellipticity vector (cross product of xivec and evec)
/// vec1: direction to first atom in path (rat[1] - rat[0])
/// vec2: direction to last atom in path (rat[nleg-1] - rat[0])
/// ndeg: path degeneracy (integer)
/// elpty: ellipticity
/// ilinit: linit + 1
double xstar(const double eps1[3], const double eps2[3],
             const double vec1[3], const double vec2[3],
             int ndeg, double elpty, int ilinit);

/// Cosine of angle between two 3-vectors.
double xxcos(const double veca[3], const double vecb[3]);

/// Helper for xstar: computes orientation-dependent factor.
double ystar(int lfin, double x, double y, double z, int iav);

} // namespace feff::genfmt
