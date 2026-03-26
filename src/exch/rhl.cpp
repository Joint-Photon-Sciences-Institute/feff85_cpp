#include "rhl.hpp"
#include "imhl.hpp"
#include <feff/constants.hpp>
#include <cmath>

namespace feff::exch {

// Interpolation coefficients for HL self-energy.
// lrs=4 panels, nrs=3 order of rs expansion, nleft=4, nright=2

// Right-branch coefficients: rcfr[lrs][nrs][nright] = [4][3][2]
// Fortran layout: rcfr(lrs, nrs, nright) with column-major
// Data listed column by column in Fortran DATA statement.
// Reordered to C row-major: rcfr[mrs][irs][j] where
//   mrs = panel index (0..3), irs = rs expansion index (0..2), j = right coeff (0..1)
static constexpr double rcfr[4][3][2] = {
    // mrs=0 (rs < 0.2)
    // rcfr[mrs][irs][j] = Fortran rcfr(mrs+1, irs+1, j+1)
    {{ -0.173963e+00, -0.116431e+00 },
     { -0.838843e-01,  0.791051e-01 },
     { -0.645803e-01, -0.628162e-01 }},
    // mrs=1 (0.2 <= rs < 1.0)
    {{ -0.173678e+00, -0.909300e-01 },
     { -0.807046e-01, -0.359401e-01 },
     { -0.731172e-01,  0.669257e-01 }},
    // mrs=2 (1.0 <= rs < 5.0)
    {{ -0.142040e+00, -0.886979e-01 },
     { -0.135577e+00, -0.379584e-01 },
     { -0.498823e-01,  0.667119e-01 }},
    // mrs=3 (rs >= 5.0)
    {{ -0.101030e+00, -0.702319e-01 },
     { -0.177556e+00, -0.419807e-01 },
     { -0.393108e-01,  0.648175e-01 }}
};

// Left-branch coefficients: rcfl[lrs][nrs][nleft] = [4][3][4]
// Same reordering from Fortran column-major.
static constexpr double rcfl[4][3][4] = {
    // mrs=0 (rs < 0.2)
    {{  0.590195e+02, -0.181726e+03,  0.184417e+03, -0.620411e+02 },
     { -0.291180e+03,  0.886023e+03, -0.895807e+03,  0.300946e+03 },
     {  0.363830e+03, -0.110486e+04,  0.111549e+04, -0.374494e+03 }},
    // mrs=1 (0.2 <= rs < 1.0)
    {{  0.478860e+01, -0.169709e+02,  0.180204e+02, -0.616427e+01 },
     { -0.926539e+01,  0.301808e+02, -0.318696e+02,  0.109158e+02 },
     {  0.460433e+01, -0.149086e+02,  0.156448e+02, -0.535127e+01 }},
    // mrs=2 (1.0 <= rs < 5.0)
    {{  0.812813e+00, -0.409425e+01,  0.450425e+01, -0.153874e+01 },
     { -0.858348e+00,  0.305836e+01, -0.345827e+01,  0.120028e+01 },
     {  0.173067e+00, -0.662794e+00,  0.749582e+00, -0.261260e+00 }},
    // mrs=3 (rs >= 5.0)
    {{  0.191145e+00, -0.173077e+01,  0.184349e+01, -0.609114e+00 },
     { -0.246947e+00,  0.743167e+00, -0.855367e+00,  0.290985e+00 },
     {  0.239738e-01, -0.100106e+00,  0.117680e+00, -0.405337e-01 }}
};

void rhl(double rs, double xk, double& erl, double& eim) {
    constexpr int nleft = 4;
    constexpr int nright = 2;

    double cleft[nleft];
    double cright[nright];

    // Calculate HL using interpolation coefficients
    double rkf = feff::fa / rs;
    double ef = rkf * rkf / 2.0;
    double wp = std::sqrt(3.0 / (rs * rs * rs));
    // Quick fix to remove jump at wp in rhl. (ala 08.01.95)
    // Use smooth transition between 2 curves in energy range dwp
    double dwp = wp / 3.0;

    int icusp;
    imhl(rs, xk, eim, icusp);

    // eim already has a factor of ef in it
    // eim also gives the position of the cusp

    double xx = xk / rkf;
    // Set to Fermi level if below Fermi level
    if (xx < 1.00001) {
        xx = 1.00001;
    }
    // Quick fix to remove jump at wp in rhl. (ala 08.01.95)
    double deltae = ((xx * xx - 1.0) * ef - wp - dwp) / dwp;

    // Calculate right hand side coefficients - determine panel index
    int mrs;
    if (rs < 0.2) {
        mrs = 0;
    } else if (rs < 1.0) {
        mrs = 1;
    } else if (rs < 5.0) {
        mrs = 2;
    } else {
        mrs = 3;
    }

    for (int j = 0; j < nright; ++j) {
        cright[j] = rcfr[mrs][0][j] * rs + rcfr[mrs][1][j] * rs * std::sqrt(rs)
                   + rcfr[mrs][2][j] * rs * rs;
    }
    double eee = -feff::pi * wp / (4.0 * rkf * ef);

    erl = 0.0;

    // Quick fix: use smooth transition
    if (icusp != 1 || std::abs(deltae) < 1.0) {
        for (int j = 0; j < nleft; ++j) {
            cleft[j] = rcfl[mrs][0][j] * rs + rcfl[mrs][1][j] * std::pow(rs, 1.5)
                       + rcfl[mrs][2][j] * rs * rs;
        }
        erl = cleft[0];
        for (int j = 1; j < nleft; ++j) {
            erl = erl + cleft[j] * std::pow(xx, j);
        }
    }

    if (icusp == 1 || std::abs(deltae) < 1.0) {
        // Right branch
        double erlr = eee / xx;
        for (int j = 0; j < nright; ++j) {
            erlr = erlr + cright[j] / std::pow(xx, j + 2);
        }
        if (std::abs(deltae) < 1.0) {
            double wr;
            if (deltae < 0.0) {
                wr = (1.0 + deltae) * (1.0 + deltae) / 2.0;
            } else {
                wr = 1.0 - (1.0 - deltae) * (1.0 - deltae) / 2.0;
            }
            erl = wr * erlr + (1.0 - wr) * erl;
        } else {
            erl = erlr;
        }
    }

    erl = erl * ef;
}

} // namespace feff::exch
