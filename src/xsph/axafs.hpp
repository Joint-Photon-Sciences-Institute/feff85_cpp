#pragma once
// Extract AXAFS from atomic cross-section.
// Converted from src/XSPH/axafs.f

#include <feff/types.hpp>
#include <feff/dimensions.hpp>

namespace feff::xsph {

/// Extract AXAFS from the atomic cross-section xsec.
/// Writes axafs.dat file.
///
/// @param em    Complex energy grid [nex]
/// @param emu   Edge energy (Hartrees)
/// @param xsec  Cross-section array [nex]
/// @param ne1   Number of horizontal grid points
/// @param ik0   Index of k=0 point
void axafs(const FeffComplex em[], double emu, const FeffComplex xsec[],
           int ne1, int ik0);

} // namespace feff::xsph
