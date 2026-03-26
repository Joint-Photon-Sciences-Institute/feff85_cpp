#pragma once
// Read user-defined energy grid from grid.inp.
// Converted from src/XSPH/rdgrid.f

#include <feff/types.hpp>
#include <feff/dimensions.hpp>

namespace feff::xsph {

/// Read grid definitions from grid.inp file.
///
/// @param em          Energy grid (complex, output for user-defined points)
/// @param ne          Current number of energy points (in/out)
/// @param nGrid       Number of grids defined (output)
/// @param iGridType   Type of each grid: 0=user, 1=energy, 2=k, 3=exponential (output)
/// @param gridMin     Minimum of each grid (output)
/// @param gridMax     Maximum of each grid (output)
/// @param gridStep    Step size of each grid (output)
/// @param nGridMax    Maximum number of grids allowed
/// @param nex         Maximum number of energy points
void rdgrid(FeffComplex em[], int& ne, int& nGrid,
            int iGridType[], double gridMin[], double gridMax[], double gridStep[],
            int nGridMax, int nex);

} // namespace feff::xsph
