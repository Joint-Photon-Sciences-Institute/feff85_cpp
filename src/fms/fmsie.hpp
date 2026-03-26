#pragma once

// Single energy FMS calculation.
// Converted from: fmsie.f
// Called by POT/scmt for self-consistent field calculations.
//
// fmsie performs full multiple scattering for one energy point,
// computing the Green's function trace for each angular momentum channel.

#include <complex>
#include <vector>
#include <Eigen/Dense>

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <feff/types.hpp>

#include "fms_types.hpp"

namespace feff::fms {

/// Full multiple scattering for a single energy point.
///
/// @param verbose  Print progress messages
/// @param iph0     Potential index of central atom (0 for absorber)
/// @param nph      Number of unique potentials
/// @param lipotx   Max angular momentum per potential [nphx+1]
/// @param ie       Energy index (1-based, for messages)
/// @param em       Complex energy (Hartree)
/// @param eref     Complex reference energy (Hartree)
/// @param ph       Phase shifts [lx+1][nphx+1] — double precision complex
/// @param rfms     FMS cluster radius (Bohr, single precision)
/// @param lfms     FMS mode: 0=crystal, 1=molecule, 2=molecule+force yprep
/// @param nat      Number of atoms in extended cluster
/// @param iphat    Potential index per atom [nat]
/// @param rath     Atom coordinates (Bohr, double) [nat][3]
/// @param[in,out] gtr  Green's function trace [lx+1][nphx+1]
///                      Accumulated: gtr(il, ip) += Tr(gg) * exp(2i*delta)/(2l+1)
/// @param data     FMS shared data (modified by yprep on first call)
void fmsie(bool verbose, int iph0, int nph, const int* lipotx,
           int ie, FeffComplex em, FeffComplex eref,
           const FeffComplex* ph,  // [(lx+1) * (nphx+1)]
           float rfms, int lfms, int nat,
           const int* iphat, const double* rath,
           Complexf* gtr,  // [(lx+1) * (nphx+1)]
           FMSData& data);

} // namespace feff::fms
