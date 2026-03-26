#pragma once
// Write diagnostic phase shift files.
// Converted from src/XSPH/wphase.f

#include <feff/types.hpp>
#include <feff/dimensions.hpp>

namespace feff::xsph {

/// Write phase data to phaseXX.dat and phminXX.dat files.
///
/// @param nph     Number of unique potentials
/// @param em      Complex energy grid [nex]
/// @param eref    Energy reference [nex][nspx]
/// @param lmax    Max angular momentum per potential [nphx+1]
/// @param ne      Total number of energy points
/// @param ph      Phase shifts [nex][2*ltot+1][nspx][nphx+1]
/// @param ntitle  Number of title lines
/// @param title   Title strings [ntitle]
void wphase(int nph, const FeffComplex em[], const FeffComplex* eref,
            const int lmax[], int ne, const FeffComplex* ph,
            int ntitle, const char title[][80]);

} // namespace feff::xsph
