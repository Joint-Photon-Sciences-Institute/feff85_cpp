#pragma once
// FMS self-consistency integration for szlz calculations.
// Converted from src/XSPH/fmssz.f

#include <feff/types.hpp>
#include <feff/dimensions.hpp>

namespace feff::xsph {

/// Perform FMS for a cluster in the self-consistency loop.
///
/// @param verbose  Enable log messages
/// @param iph0    Central potential index
/// @param ie      Energy point index
/// @param em      Complex energy
/// @param eref    Energy reference
/// @param ph      Phase shifts [lx+1][nphx+1]
/// @param nph     Number of unique potentials
/// @param rfms    FMS radius
/// @param lfms    FMS angular momentum flag
/// @param nat     Number of atoms
/// @param iphat   Potential type per atom [natx]
/// @param rat     Atomic coordinates [3][natx]
/// @param amat    Angular coefficient matrix
/// @param lipotx  Max l per potential [nphx+1]
/// @param gctr    Central atom Green's function (real part) [2][2][3][lx+1][nphx+1]
/// @param gtr     Full Green's function (complex) [2][2][3][lx+1][nphx+1]
void fmssz(bool verbose, int iph0, int ie, FeffComplex em, FeffComplex eref,
           const FeffComplex* ph, int nph,
           float rfms, int lfms, int nat, const int iphat[],
           const double* rat, const float amat[], const int lipotx[],
           float* gctr, FeffComplex* gtr);

} // namespace feff::xsph
