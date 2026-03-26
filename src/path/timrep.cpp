// Time-reversal and standard ordering of paths.
// Converted from: src/PATH/timrep.f

#include "timrep.hpp"
#include "mpprmp.hpp"
#include "phash.hpp"
#include <feff/dimensions.hpp>

namespace feff::path {

void timrep(int npat, int ipat[], float rx[], float ry[], float rz[],
            double& dhash, int ipol, int ispin,
            const double evec[3], const double xivec[3], int eels,
            const AtomData& atoms) {

    // Prepare icase for mpprmp
    int icase = -1;
    if (eels == 1) icase = 7;

    int nleg = npat + 1;
    ipat[npat] = 0;  // central atom at end

    // Initialize coordinate arrays
    for (int i = 0; i < npatx; ++i) {
        rx[i] = 0.0f;
        ry[i] = 0.0f;
        rz[i] = 0.0f;
    }

    // Compute standard-frame coordinates and hash for this ordering
    mpprmp(npat, ipat, rx, ry, rz, ipol, ispin, evec, xivec, icase, atoms);
    phash(npat, ipat, rx, ry, rz, dhash, atoms);

    // If single-bounce (npat <= 1), path is symmetrical -- done
    if (npat <= 1) return;

    // Make time-reversed path
    int ipat0[npatx + 1];
    float rx0[npatx], ry0[npatx], rz0[npatx];

    ipat0[nleg - 1] = ipat[nleg - 1];  // central atom
    for (int i = 0; i < npat; ++i) {
        ipat0[i] = ipat[nleg - 2 - i];  // reverse order
    }
    for (int i = 0; i < npatx; ++i) {
        rx0[i] = 0.0f;
        ry0[i] = 0.0f;
        rz0[i] = 0.0f;
    }

    mpprmp(npat, ipat0, rx0, ry0, rz0, ipol, ispin, evec, xivec, icase, atoms);
    double dhash0;
    phash(npat, ipat0, rx0, ry0, rz0, dhash0, atoms);

    // Turn off path reversal for spin-polarized cases
    if (ispin != 0 && ipol != 0) dhash0 = dhash + 1.0;

    // Want representation with smallest hash number
    if (dhash0 < dhash) {
        dhash = dhash0;
        for (int i = 0; i < npat; ++i) {
            ipat[i] = ipat0[i];
            rx[i] = rx0[i];
            ry[i] = ry0[i];
            rz[i] = rz0[i];
        }
    }
}

} // namespace feff::path
