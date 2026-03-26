#pragma once
// Nuclear potential construction and radial mesh setup.
// Converted from: nucdev.f

#include "atom_types.hpp"

namespace feff::atom {

/// Construct nuclear potential and set up radial mesh.
///
/// av[max_dev]  — development coefficients at origin of nuclear potential (output)
/// dr[atom_grid] — radial mesh points (output)
/// dv[atom_grid] — nuclear potential on mesh (output)
/// dz           — nuclear charge
/// hx           — exponential mesh step
/// nuc          — nuclear radius index (1 = point charge) (in/out)
/// np           — number of tabulation points
/// ndor         — number of development coefficients (must be >= 5)
/// dr1          — first mesh point (in/out, may be adjusted)
///
/// All arrays are 0-based in C++.
void nucdev(double av[], double dr[], double dv[], double dz, double hx,
            int& nuc, int np, int ndor, double& dr1);

} // namespace feff::atom
