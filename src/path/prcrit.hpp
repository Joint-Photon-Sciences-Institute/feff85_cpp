#pragma once
// Prepare scattering amplitude arrays for pathfinder criteria.
// Converted from: src/PATH/prcrit.f

#include "path_data.hpp"
#include <string>

namespace feff::path {

/// Prepare fbeta, ckspc, etc. arrays for path criteria.
/// Reads phase.pad via feff::common::read_xsph().
/// neout: number of energy grid points (output).
/// ik0out: index of k=0 point (output).
/// cksp[0..ne-1]: |p| at each energy point (single precision).
/// fbeta[fbeta_dim]: |f(beta)| for each angle, pot, energy (single precision).
/// ckspc[0..nncrit-1]: |p| at criterion points.
/// fbetac[fbetac_dim]: |f(beta)| at criterion points.
/// potlb0[0..nphx]: potential labels.
/// xlam[0..ne-1]: mean free path in Angstrom.
/// xlamc[0..nncrit-1]: mfp at criterion points.
/// nncrit: number of criterion points used (output).
void prcrit(int& neout, int& nncrit, int& ik0out,
            float cksp[], float fbeta[],
            float ckspc[], float fbetac[],
            std::string potlb0[],
            float xlam[], float xlamc[]);

} // namespace feff::path
