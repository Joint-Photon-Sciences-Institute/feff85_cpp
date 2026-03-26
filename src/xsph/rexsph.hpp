#pragma once
// Read XSPH module configuration from JSON.
// Converted from src/XSPH/rexsph.f

#include <feff/types.hpp>
#include <feff/dimensions.hpp>
#include <string>

namespace feff::xsph {

/// Module configuration parameters read from xsph.json.
struct XsphConfig {
    int mphase = 0;
    int ipr2 = 0;
    int ixc = 0;
    int ixc0 = 0;
    int ispec = 0;
    int lreal = 0;
    int lfms2 = 0;
    int l2lp = 0;
    int iPl = 0;
    int iGrid = 0;
    int izstd = 0;
    int ifxc = 0;
    int ipmbse = 0;
    int itdlda = 0;
    int nonlocal = 0;
    int ibasis = 0;

    double vr0 = 0.0;
    double vi0 = 0.0;
    double rgrd = 0.05;
    float rfms2 = 0.0f;
    double gamach = 0.0;
    double xkstep = 0.0;
    double xkmax = 0.0;
    double vixan = 0.0;

    int nph = 0;
    int nat = 0;
    int ipol = 0;
    int ispin = 0;
    int le2 = 0;
    double angks = 0.0;
    FeffComplex ptz[3][3] = {};

    int lmaxph[nphx + 1] = {};
    double spinph[nphx + 1] = {};
    char potlbl[nphx + 1][7] = {};
    int iatph[nphx + 1] = {};
    int iphat[natx] = {};
    double rat[3][natx] = {};
};

/// Read xsph module configuration from JSON files.
///
/// @param config  Output configuration structure
void rexsph(XsphConfig& config);

} // namespace feff::xsph
