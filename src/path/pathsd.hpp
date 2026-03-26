#pragma once
// Path degeneracy checker: eliminates degenerate paths and writes paths.dat.
// Converted from: src/PATH/pathsd.f

#include "path_data.hpp"
#include <string>

namespace feff::path {

/// Eliminate path degeneracies, compute importance factors, write paths.dat.
/// Reads paths.bin (written by paths()), writes paths.dat and crit.dat.
void pathsd(const float ckspc[], const float fbetac[], const float xlamc[],
            int ne, int ik0, const float cksp[],
            const float fbeta[], const float xlam[],
            float critpw, int ipr2, int nncrit,
            const std::string potlbl[],
            int ipol, int ispin, const double evec[3], const double xivec[3],
            int eels);

} // namespace feff::path
