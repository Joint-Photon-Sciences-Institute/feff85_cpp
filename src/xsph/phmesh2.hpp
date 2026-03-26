#pragma once
// Alternative energy mesh with user-defined grid support.
// Converted from src/XSPH/phmesh2.f
//
// Supports both FEFF84-style grids and user-defined grids
// read from grid.inp.

#include <feff/types.hpp>
#include <feff/dimensions.hpp>

namespace feff::xsph {

/// Generate energy mesh with user-defined grid support.
/// See phmesh.hpp for parameter descriptions (identical interface).
/// @param iGrid  0 = use FEFF84 grids, nonzero = read from grid.inp
void phmesh2(int iprint, int ispec, double edge, double emu,
             double vi0, double gamach,
             double xkmax, double xkstep, double vixan,
             int& ne, int& ne1, FeffComplex em[], int& ik0, int& ne3,
             int iGrid);

} // namespace feff::xsph
