// Grid interpolation — converted from src/COMMON/fixdsp.f, fixdsx.f, fixvar.f

#include "grid_interpolation.hpp"
#include "logging.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <cmath>
#include <algorithm>

// Forward declaration of terp from math module
namespace feff::math {
    void terp(const double x[], const double y[], int n, int m,
              double x0, double& y0);
}

namespace feff::common {

// Local grid helpers (1-based j, like Fortran)
static constexpr double xx00 = 8.8;

static inline double xxx(int j, double delta) {
    return -xx00 + (j - 1) * delta;
}

static inline double rrr(int j, double delta) {
    return std::exp(-xx00 + (j - 1) * delta);
}

static inline int jjj(double r, double delta) {
    return static_cast<int>((std::log(r) + xx00) / delta) + 1;
}

void fixdsp(double dxorg, double dxnew,
            const double dgc0[251], const double dpc0[251],
            double dgcx[], double dpcx[], int& jnew) {

    // Find last non-zero point in original arrays (1-based)
    int imax = 0;
    for (int i = 250; i >= 0; --i) {  // 0-based: 250..0 → Fortran 251..1
        if (std::abs(dgc0[i]) >= 1.0e-11 || std::abs(dpc0[i]) >= 1.0e-11) {
            imax = i + 1;  // convert to 1-based
            break;
        }
    }
    if (imax == 0) {
        logger().wlog(" Should never see this line from sub fixdsp");
    }

    int jmax = std::min(imax, 250) + 1;  // 1-based

    // Build original grid x-values
    double xorg[nrptx];
    double delta = dxorg;
    for (int j = 1; j <= jmax; ++j) {
        xorg[j - 1] = xxx(j, delta);
    }
    double rmax = rrr(jmax, delta);

    // Build new grid x-values
    double xnew[nrptx];
    delta = dxnew;
    jnew = jjj(rmax, delta);
    for (int j = 1; j <= jnew; ++j) {
        xnew[j - 1] = xxx(j, delta);
    }

    // Interpolate to new grid using cubic interpolation (m=3)
    for (int j = 1; j <= jnew; ++j) {
        feff::math::terp(xorg, dgc0, jmax, 3, xnew[j - 1], dgcx[j - 1]);
        feff::math::terp(xorg, dpc0, jmax, 3, xnew[j - 1], dpcx[j - 1]);
    }

    // Zero arrays past rmax
    for (int j = jnew; j < nrptx; ++j) {
        dgcx[j] = 0.0;
        dpcx[j] = 0.0;
    }
}

void fixdsx(int iph, double dxorg, double dxnew,
            const double* dgc, const double* dpc,
            double* dgcn, double* dpcn) {
    // fixdsx interpolates 30 orbitals for one potential type.
    // The Fortran arrays are dgc(251, 30, 0:nphx+1) and dpc similarly.
    // We treat dgc/dpc as flat arrays and compute offsets.

    // Strides for dgc/dpc: dgc(i, iorb, iph) → dgc[(iph)*(251*30) + (iorb)*251 + i]
    constexpr int dim1 = 251;
    constexpr int dim2 = 30;
    constexpr int stride_orb = dim1;
    constexpr int stride_ph = dim1 * dim2;

    // Strides for dgcn/dpcn: dgcn(j, iorb) → dgcn[iorb * nrptx + j]
    constexpr int nstride = nrptx;

    for (int iorb = 0; iorb < 30; ++iorb) {
        // Extract single orbital from the 3D array
        double dgc0[251], dpc0[251];
        bool all_zero = true;
        for (int i = 0; i < 251; ++i) {
            dgc0[i] = dgc[iph * stride_ph + iorb * stride_orb + i];
            dpc0[i] = dpc[iph * stride_ph + iorb * stride_orb + i];
            if (std::abs(dgc0[i]) >= 1.0e-11 || std::abs(dpc0[i]) >= 1.0e-11)
                all_zero = false;
        }

        double dgcx[nrptx], dpcx[nrptx];
        int jnew;
        if (all_zero) {
            // Empty orbital slot - zero the output and set jnew = 0
            jnew = 0;
            for (int j = 0; j < nrptx; ++j) {
                dgcx[j] = 0.0;
                dpcx[j] = 0.0;
            }
        } else {
            fixdsp(dxorg, dxnew, dgc0, dpc0, dgcx, dpcx, jnew);
        }

        // Copy results into output arrays
        for (int j = 0; j < nrptx; ++j) {
            dgcn[iorb * nstride + j] = dgcx[j];
            dpcn[iorb * nstride + j] = dpcx[j];
        }
    }
}

void fixvar(double rmt, const double edens[251], const double vtot[251],
            const double dmag[251], double vint, double rhoint,
            double dxorg, double dxnew, int jumprm,
            double& vjump, double ri[], double vtotph[],
            double rhoph[], double dmagx[]) {

    double xorg[nrptx], xnew[nrptx];

    // Build original grid
    double delta = dxorg;
    int jmtorg = jjj(rmt, delta);
    int jriorg = jmtorg + 1;
    int jrior1 = jriorg + 1;
    for (int j = 1; j <= jrior1; ++j) {
        xorg[j - 1] = xxx(j, delta);
    }

    // Build new grid
    delta = dxnew;
    int jmtnew = jjj(rmt, delta);
    int jrinew = jmtnew + 1;
    int jrine1 = jrinew + 1;
    for (int j = 1; j <= jrine1; ++j) {
        xnew[j - 1] = xxx(j, delta);
    }

    // Interpolate inside muffin tin
    for (int j = 1; j <= jrinew; ++j) {
        feff::math::terp(xorg, vtot, jriorg, 3, xnew[j - 1], vtotph[j - 1]);
        feff::math::terp(xorg, edens, jrior1, 3, xnew[j - 1], rhoph[j - 1]);
        feff::math::terp(xorg, dmag, jrior1, 3, xnew[j - 1], dmagx[j - 1]);
    }

    // Handle potential jump at MT boundary
    if (jumprm == 1) {
        double xmt = std::log(rmt);
        double vmt;
        feff::math::terp(xorg, vtot, jriorg, 3, xmt, vmt);
        vjump = vint - vmt;
    }
    if (jumprm > 0) {
        for (int j = 0; j < jrinew; ++j) {
            vtotph[j] += vjump;
        }
    }

    // Set radial grid and fill interstitial region
    delta = dxnew;
    for (int j = 1; j <= nrptx; ++j) {
        ri[j - 1] = rrr(j, delta);
    }
    for (int j = 0; j < jrinew; ++j) {
        rhoph[j] /= (4.0 * pi);
    }
    for (int j = jrinew; j < nrptx; ++j) {
        vtotph[j] = vint;
        rhoph[j] = rhoint / (4.0 * pi);
        dmagx[j] = 0.0;
    }
}

} // namespace feff::common
