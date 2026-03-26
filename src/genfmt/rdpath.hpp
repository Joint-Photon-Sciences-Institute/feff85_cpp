#pragma once

// Read path data from paths.dat file.
// Converted from GENFMT/rdpath.f

#include <feff/dimensions.hpp>
#include <string>
#include <iostream>

namespace feff::genfmt {

/// Read a single path from the paths.dat file.
///
/// in: input stream (already open)
/// done: set to true when no more paths remain
/// ipol: polarization flag
/// potlbl: potential labels (updated for each atom)
/// rat: atom positions [3][legtot+2] (output, converted from Angstrom to Bohr)
/// ri: leg lengths (output)
/// beta: Euler beta angles (output)
/// eta: combined Euler angles (output)
/// deg: path degeneracy (output)
/// ipot: potential indices (output)
/// nsc: number of scatterers (output, = nleg - 1)
/// nleg: number of legs (output)
/// npot: number of potentials (input, for range checking)
/// ipath: path index (output)
void rdpath(std::istream& in, bool& done, int ipol,
            std::string potlbl[], double rat[][legtot + 2],
            double ri[], double beta[], double eta[],
            double& deg, int ipot[], int& nsc, int& nleg,
            int npot, int& ipath);

} // namespace feff::genfmt
