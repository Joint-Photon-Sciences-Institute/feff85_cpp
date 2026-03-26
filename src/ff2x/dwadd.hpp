#pragma once
// Debye-Waller factor addition to paths.
// Converted from: src/FF2X/ff2gen.f (dwadd subroutine)
// Adds DW factors, cumulants, and sums path contributions.

#include "ff2x_types.hpp"

#include <complex>
#include <vector>

namespace feff::ff2x {

/// Add Debye-Waller factors and sum path contributions into cchi.
/// Replaces Fortran subroutine dwadd.
///
/// @param ntotal   Number of paths in path_list
/// @param nptot    Total paths in feff.pad
/// @param p        FF2X parameters
/// @param path_list  List of paths with user sigma^2
/// @param pad      Data from feff.pad
/// @param xs       Data from xsect.bin
/// @param nkx      Number of fine k-grid points
/// @param xk0      Fine k-grid (interpolation)
/// @param xkp      Fine k-grid (output)
/// @param cchi     Accumulated complex chi (nkx points, modified in place)
/// @param iabs     Current absorber index
/// @param nused    Output: number of paths actually used
void dwadd(int ntotal, int nptot, const FF2xParams& p,
           const std::vector<PathListEntry>& path_list,
           FeffPadData& pad, const XsectData& xs,
           int nkx, const double xk0[], const double xkp[],
           FeffComplex cchi[], int iabs, int& nused);

} // namespace feff::ff2x
