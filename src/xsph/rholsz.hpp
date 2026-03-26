#pragma once
// Spectroscopic density of states (s_z version, no projection).
// Converted from src/XSPH/rholsz.f

#include <feff/types.hpp>
#include <feff/dimensions.hpp>

namespace feff::xsph {

/// Calculate spectroscopic DOS without orbital projection.
/// Similar to rholat but simpler (no icount, no iorb projection).
///
/// See rholat.hpp for parameter descriptions.
void rholsz(double dx, double x0, const double ri[],
            FeffComplex em, int ixc, double rmt, double rnrm,
            const double vtot[], const double vvalgs[],
            const double xnval[],
            const double dgcn[][30], const double dpcn[][30],
            FeffComplex eref,
            const double adgc[][30], const double adpc[][30],
            FeffComplex* xrhole, FeffComplex* xrhoce,
            FeffComplex ph[], int iz, double xion, int iunf,
            int ihole, int lmaxsc);

} // namespace feff::xsph
