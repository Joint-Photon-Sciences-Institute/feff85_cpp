#pragma once
// Read input for the PATH module (geometry, global, path parameters).
// Converted from: src/PATH/repath.f

#include <string>

namespace feff::path {

/// Read all input needed by the pathfinder.
/// Reads geom.json, global.json, and path.json.
/// All outputs are set from JSON input files.
void repath(int& ms, int& mpath, int& ipr4,
            float& pcritk, float& pcrith, int& nncrit,
            float& rmax, int& nlegxx, float& rfms2, float& critpw,
            int& nat, double rat[][3], int iphat[], int ibounc[],
            int& ipol, int& ispin, double evec[3], double xivec[3]);

} // namespace feff::path
