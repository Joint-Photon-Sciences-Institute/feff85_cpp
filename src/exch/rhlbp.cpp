#include "rhlbp.hpp"
#include "quinn.hpp"
#include <feff/constants.hpp>
#include "../math/interpolation.hpp"
#include <cmath>
#include <fstream>
#include <stdexcept>

namespace feff::exch {

// 2D linear interpolation and extrapolation (analog of terp.f in 2D).
// Internal helper — not exposed in header.
static void terp2d(const double x[], const double y[], const double z[],
                   int nx, int ny, double x0, double y0, double& z0) {
    // Find out between which x points x0 lies
    int ix = feff::math::locat(x0, nx, x);
    if (ix < 1) ix = 1;
    if (ix > nx - 1) ix = nx - 1;
    // Convert to 0-based for C++ array access
    int ix0 = ix - 1;  // locat returns 1-based index
    if (x[ix0 + 1] - x[ix0] == 0.0) {
        throw std::runtime_error("TERP-1: zero x interval in terp2d");
    }

    // Find out between which y points y0 lies
    int iy = feff::math::locat(y0, ny, y);
    if (iy < 1) iy = 1;
    if (iy > ny - 1) iy = ny - 1;
    int iy0 = iy - 1;
    if (y[iy0 + 1] - y[iy0] == 0.0) {
        throw std::runtime_error("TERP-1: zero y interval in terp2d");
    }

    double dx = (x0 - x[ix0]) / (x[ix0 + 1] - x[ix0]);
    double dy = (y0 - y[iy0]) / (y[iy0 + 1] - y[iy0]);

    // z is stored as z[ix * ny + iy] (row-major, Fortran was z(nx, ny))
    // Fortran: z(ix, iy) => C: z[(ix-1)*ny + (iy-1)]
    // Note: reproducing exact Fortran behavior where z1 == z2 (apparent bug in original)
    double z1 = z[ix0 * ny + iy0] + dx * (z[(ix0 + 1) * ny + iy0] - z[ix0 * ny + iy0]);
    double z2 = z[ix0 * ny + iy0] + dx * (z[(ix0 + 1) * ny + iy0] - z[ix0 * ny + iy0]);
    z0 = z1 + dy * (z2 - z1);
}

// Dimensions for bphl.dat grid
static constexpr int nrs_bp = 21;
static constexpr int nx_bp = 51;

// Cached data from bphl.dat
static bool initialized = false;
static double rsmesh[nrs_bp];
static double xmesh[nx_bp];
static double sigma_re[nrs_bp * nx_bp];  // sigma(irs, ik, 1) -> real part
static double sigma_im[nrs_bp * nx_bp];  // sigma(irs, ik, 2) -> imaginary part

void rhlbp(double rs, double xk, double& erl, double& eim) {
    double xf = feff::fa / rs;
    double ef = xf * xf / 2.0;
    double wp = std::sqrt(3.0 / (rs * rs * rs)) / ef;
    double xk0 = xk / xf;
    double xx = (xk0 * xk0 - 1.0) / std::sqrt(rs);

    if (!initialized) {
        // Read self energy for grid points from bphl.dat
        std::ifstream infile("bphl.dat");
        if (!infile.is_open()) {
            throw std::runtime_error("rhlbp: cannot open bphl.dat");
        }

        xmesh[0] = 0.0;
        for (int irs = 0; irs < nrs_bp; ++irs) {
            // sigma(irs, 0, re/im) = 0.0
            sigma_re[irs * nx_bp + 0] = 0.0;
            sigma_im[irs * nx_bp + 0] = 0.0;
            for (int ik = 1; ik < nx_bp; ++ik) {
                double rs_read, x_read, sig_re, sig_im;
                infile >> rs_read >> x_read >> sig_re >> sig_im;
                rsmesh[irs] = rs_read;
                xmesh[ik] = x_read;
                sigma_re[irs * nx_bp + ik] = sig_re;
                sigma_im[irs * nx_bp + ik] = sig_im;
            }
        }
        initialized = true;
        infile.close();
    }

    terp2d(rsmesh, xmesh, sigma_re, nrs_bp, nx_bp, rs, xx, erl);
    terp2d(rsmesh, xmesh, sigma_im, nrs_bp, nx_bp, rs, xx, eim);

    // Transfer to atomic units
    erl = erl / rs / feff::hart;
    eim = eim / rs / feff::hart;

    double ei;
    quinn(xk0, rs, wp, ef, ei);
    if (eim >= ei) eim = ei;
}

} // namespace feff::exch
