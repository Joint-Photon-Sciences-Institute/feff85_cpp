#pragma once
// Build polarization tensor and rotate coordinates.
// Converted from: src/RDINP/mkptz.f

#include <feff/types.hpp>

namespace feff::rdinp {

/// Build the polarization tensor ptz and optionally rotate atom coordinates
/// so that the z-axis aligns with the spin or propagation direction.
///
/// @param ipol   Polarization type: 0=random, 1=linear/elliptical, 2=circular
/// @param elpty  Ellipticity parameter (used when ipol=1)
/// @param evec   Polarization vector (3), may be modified
/// @param xivec  Incidence/propagation direction (3), may be modified
/// @param ispin  Spin flag
/// @param spvec  Spin vector (3), may be modified
/// @param nat    Number of atoms
/// @param rat    Atom coordinates [nat][3], rotated in place
/// @param angks  Output: angle between k-vector and spin
/// @param le2    In/out: 0=E1 only, 1=E1+M1, 2=E1+E2, 3=E1+E2+M1
/// @param ptz    Output: polarization tensor [3][3], indexed as [i+1][j+1]
///
/// Replaces Fortran subroutine mkptz.
void mkptz(int ipol, double elpty, double evec[3], double xivec[3],
           int ispin, double spvec[3], int nat, double rat[][3],
           double& angks, int& le2, FeffComplex ptz[3][3]);

/// Rotate a 3-vector by Euler angles defined by cos/sin of theta and phi.
/// Replaces Fortran subroutine rotate(vec, cst, snt, csf, snf).
void rotate(double vec[3], double cst, double snt, double csf, double snf);

} // namespace feff::rdinp
