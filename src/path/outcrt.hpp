#pragma once
// Output criterion: plane-wave importance factor and pathfinder criteria.
// Converted from: src/PATH/outcrt.f

#include "path_data.hpp"

namespace feff::path {

/// Compute plane-wave importance factor and recalculate pathfinder criteria.
/// xport: integrated pw importance factor.
/// xheap: heap criterion for this path.
/// xheapr: heap criterion for time-reversed path.
/// xout: keep criterion.
/// xcalcx: running max criterion (passed in/out).
void outcrt(int npat, const int ipat[],
            const float ckspc[], int nncrit,
            const float fbetac[], const float xlamc[],
            int ne, int ik0, const float cksp[],
            const float fbeta[], const float xlam[],
            const int ipot[],
            float& xport, float& xheap, float& xheapr,
            float& xout, float& xcalcx,
            const AtomData& atoms);

} // namespace feff::path
