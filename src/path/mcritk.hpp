#pragma once
// Keep path criterion: decides whether to output a path.
// Converted from: src/PATH/mcritk.f

#include "path_data.hpp"

namespace feff::path {

/// Compute output importance factor xout for a complete path.
/// Returns xout = -1 if last atom is central atom (undefined).
/// xcalcx: max xcalc encountered so far. Set to -1 to reset;
///         computed only on first call (NN SS path).
void mcritk(int npat, const int ipat[], const float ri[], const float beta[],
            const int indbet[], const int ipot[], int nncrit,
            const float fbetac[], const float xlamc[], const float ckspc[],
            float& xout, float& xcalcx);

} // namespace feff::path
