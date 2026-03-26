// Set radial grid indices for MT and Norman radii.
// Converted from src/POT/sidx.f

#include "sidx.hpp"
#include "../common/radial_grid.hpp"
#include "../common/logging.hpp"
#include <sstream>
#include <iomanip>

namespace feff::pot {

void sidx(const double* rholap, int npts, double& rmt, double& rnrm,
          int& imax, int& imt, int& inrm)
{
    using feff::common::ii;
    using feff::common::rr;

    imt = ii(rmt);
    inrm = ii(rnrm);

    // Set imax (last non-zero rholap data)
    // rholap is 0-based; Fortran loop from imt to npts (1-based)
    imax = imt;
    for (int i = imt; i <= npts; ++i) {
        // i is 1-based grid index; rholap[i-1] is the 0-based array element
        if (rholap[i - 1] <= 1.0e-5) break;
        imax = i;
    }

    // Move Norman radius if density is zero inside rnrm
    if (inrm > imax) {
        inrm = imax;
        rnrm = rr(inrm);
        std::ostringstream slog;
        slog << " Moved rnrm.  New rnrm (au) " << std::scientific
             << std::setprecision(5) << rnrm;
        feff::common::logger().wlog(slog.str());
    }
    if (imt > imax) {
        imt = imax;
        rmt = rr(imt);
        std::ostringstream slog;
        slog << " Moved rmt.  New rmt (au) " << std::scientific
             << std::setprecision(5) << rmt;
        feff::common::logger().wlog(slog.str());
    }
}

} // namespace feff::pot
