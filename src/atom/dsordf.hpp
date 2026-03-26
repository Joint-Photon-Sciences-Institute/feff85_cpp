#pragma once
// Differential integral calculation using Simpson's method.
// Converted from: src/ATOM/dsordf.f
//
// Computes various overlap and norm integrals of radial wave functions.
// The workspace arrays dg/ag/dp/ap in DiracWorkspaceReal are used for
// data exchange with calling routines (ortdat, etc.).

#include "atom_types.hpp"

namespace feff::atom {

/// Calculate differential integrals by Simpson's method.
///
/// Integrates hg * r^n where hg is constructed from orbital wave functions
/// according to the jnd flag:
///   jnd=1:  hg(l) = cg(l,i)*cg(l,j) + cp(l,i)*cp(l,j)
///   jnd=-1: same as jnd=1, multiplied by dg
///   jnd=2:  hg(l) = cg(l,i)*cp(l,j)
///   jnd=-2: same as jnd=2, multiplied by dg
///   jnd=3:  hg(l) = dg(l)*cg(l,i) + dp(l)*cp(l,j)
///   jnd=4:  hg(l) = dg(l)*dg(l) + dp(l)*dp(l)
///   jnd>=5: hg constructed by calling program
///
/// @param i    First orbital index (1-based)
/// @param j    Second orbital index or max point (1-based)
/// @param n    Power of r in the integrand
/// @param jnd  Type of integral (see above)
/// @param a    Power behavior at origin for dg/dp/hg
/// @param state Aggregate atom state
/// @return     Value of the integral
double dsordf(int i, int j, int n, int jnd, double a, AtomState& state);

} // namespace feff::atom
