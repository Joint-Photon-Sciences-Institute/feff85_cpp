// Convert experimental units (Angstrom, eV) to atomic units (bohr, hartree).
// Converted from: src/RDINP/eu2au.f

#include "eu2au.hpp"
#include <feff/constants.hpp>

namespace feff::rdinp {

void eu2au(FeffInput& inp, int nat, double rat[][3]) {
    // Distance conversions (Angstrom -> bohr)
    inp.rmax  = static_cast<float>(inp.rmax / bohr);
    inp.rfms1 = static_cast<float>(inp.rfms1 / bohr);
    inp.rfms2 = inp.rfms2 / bohr;
    inp.rdirec = static_cast<float>(inp.rdirec / bohr);

    // Energy conversions (eV -> hartree)
    inp.vr0    = inp.vr0 / hart;
    inp.vi0    = inp.vi0 / hart;
    inp.vrcorr = inp.vrcorr / hart;
    inp.vicorr = inp.vicorr / hart;
    inp.gamach = inp.gamach / hart;
    inp.ecv    = inp.ecv / hart;
    inp.emin   = inp.emin / hart;
    inp.emax   = inp.emax / hart;
    inp.eimag  = inp.eimag / hart;
    inp.vixan  = inp.vixan / hart;

    // k-space conversions (1/Angstrom -> 1/bohr)
    inp.xkstep = inp.xkstep * bohr;
    inp.xkmax  = inp.xkmax * bohr;

    // Volume conversion (Angstrom^3 -> bohr^3)
    inp.totvol = inp.totvol / (bohr * bohr * bohr);

    // Convert atom positions (Angstrom -> bohr)
    // Note: Fortran uses 1-based indexing; C++ uses 0-based
    for (int iat = 0; iat < nat; ++iat) {
        for (int i = 0; i < 3; ++i) {
            rat[iat][i] = rat[iat][i] / bohr;
        }
    }

    // Convert overlap radii (Angstrom -> bohr)
    for (int iph = 0; iph <= inp.nph; ++iph) {
        for (int iovr = 0; iovr < inp.novr[iph]; ++iovr) {
            inp.rovr[iph][iovr] = inp.rovr[iph][iovr] / bohr;
        }
    }
}

} // namespace feff::rdinp
