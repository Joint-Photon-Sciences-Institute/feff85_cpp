#pragma once
// Main pathfinder: finds multiple scattering paths using heap algorithm.
// Converted from: src/PATH/paths.f

#include "path_data.hpp"

namespace feff::path {

/// Find multiple scattering paths.
/// All coordinates are in single precision Angstroms.
/// pcrith: cut-off fraction for heap (building paths).
/// pcritk: cut-off fraction for output (keeping paths).
/// critpw: plane-wave criterion percentage.
/// nncrit: number of criterion energy points.
/// rmax: max half-path length (one-way distance, doubled internally).
/// nlegxx: max number of legs user wants to consider.
/// rfms: FMS cluster radius.
/// nat: number of atoms (1-based Fortran count; converted to 0-based internally).
/// ratdp[1..nat][3]: atom positions in double precision.
/// iphat[1..nat]: potential type for each atom.
/// ibounc[1..nat]: first-bounce flag/degeneracy for each atom.
void paths(const float ckspc[], const float fbetac[], const float xlamc[],
           float pcritk, float pcrith, float critpw,
           int nncrit, float& rmax, int nlegxx, float rfms,
           int& nat, const double ratdp[][3], const int iphat[], const int ibounc[]);

} // namespace feff::path
