// Read xsect.bin (cross-section data).
// Converted from: src/FF2X/ff2gen.f (rdxbin subroutine)
// Uses read_xsect_json instead of the original binary/text format.

#include "rdxbin.hpp"

#include <feff/constants.hpp>
#include <feff/json_io.hpp>

#include <cmath>

namespace feff::ff2x {

// Compute signed k from energy (same as Fortran getxk)
static double getxk(double e) {
    if (e >= 0.0) return std::sqrt(2.0 * e);
    return -std::sqrt(-2.0 * e);
}

XsectData read_xsect_bin(double s02_in, int mbconv) {
    XsectData xs;

    // Read from xsect.json (replaces the old xsect.bin text read)
    double er[nex]{}, ei[nex]{}, xsn[nex]{}, col4[nex]{}, col5[nex]{};
    double s02p_raw = 0.0, erelax_raw = 0.0, wp_raw = 0.0, edgep_raw = 0.0, emu_raw = 0.0;
    double gamach_raw = 0.0;
    int ne_raw = 0, ne1_raw = 0, ik0_raw = 0;
    int ntit_raw = 0;
    std::string titles_raw[nheadx];

    json_io::read_xsect_json(ntit_raw, titles_raw,
                             s02p_raw, erelax_raw, wp_raw,
                             edgep_raw, emu_raw, gamach_raw,
                             ne_raw, ne1_raw, ik0_raw,
                             er, ei, xsn, col4, col5);

    xs.s02p   = s02p_raw;
    xs.erelax = erelax_raw;
    xs.wp     = wp_raw;
    xs.edgep  = edgep_raw;
    xs.gamach  = gamach_raw / hart;  // Convert eV to code units
    xs.ne1    = ne1_raw;
    xs.ik0    = ik0_raw;
    xs.nxsec  = ne_raw;
    xs.ntitle = ntit_raw;

    // Override s02 if mbconv > 0 or s02 too small (matching Fortran rdxbin)
    double s02_use = s02_in;
    if (mbconv > 0 || s02_in <= 0.1) s02_use = s02p_raw;
    xs.s02_eff = s02_use;
    // Note: s02_use is returned via the caller's parameter; here we
    // just store the raw value. The caller handles the override.

    for (int i = 0; i < ntit_raw; ++i) {
        xs.title[i] = titles_raw[i];
    }

    // Fill arrays, converting to code units
    for (int i = 0; i < ne_raw; ++i) {
        xs.xsec[i] = FeffComplex(col4[i], col5[i]);
        xs.emxs[i] = FeffComplex(er[i], ei[i]) / hart;  // eV to Hartree
        xs.xkxs[i] = getxk(xs.emxs[i].real() - xs.edgep);
        xs.omega[i] = xs.emxs[i].real() - xs.edgep + emu_raw;
        xs.xsnorm[i] = xsn[i];
    }

    xs.mbconv = mbconv;

    return xs;
}

} // namespace feff::ff2x
