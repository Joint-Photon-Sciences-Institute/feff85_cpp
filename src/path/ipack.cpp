// Integer packing/unpacking for path atom indices.
// Converted from: src/PATH/ipack.f

#include "ipack.hpp"
#include "../par/parallel.hpp"

namespace feff::path {

static constexpr int ifac1 = 1290;
static constexpr int ifac2 = 1290 * 1290;

void ipack(int iout[3], int n, const int ipat[]) {
    if (n > 8) feff::par::par_stop("ipack n too big");

    int itmp[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    for (int i = 0; i < n; ++i) {
        itmp[i] = ipat[i];
    }

    iout[0] = n       + itmp[0] * ifac1 + itmp[1] * ifac2;
    iout[1] = itmp[2] + itmp[3] * ifac1 + itmp[4] * ifac2;
    iout[2] = itmp[5] + itmp[6] * ifac1 + itmp[7] * ifac2;
}

void upack(const int iout[3], int& n, int ipat[]) {
    int nmax = n;
    if (nmax > 8) feff::par::par_stop("nmax > 8 in upack");

    n = iout[0] % ifac1;
    if (n > nmax) feff::par::par_stop("nmax in upack too small");

    int itmp[8];
    itmp[0] = (iout[0] % ifac2) / ifac1;
    itmp[1] = iout[0] / ifac2;
    itmp[2] = iout[1] % ifac1;
    itmp[3] = (iout[1] % ifac2) / ifac1;
    itmp[4] = iout[1] / ifac2;
    itmp[5] = iout[2] % ifac1;
    itmp[6] = (iout[2] % ifac2) / ifac1;
    itmp[7] = iout[2] / ifac2;

    for (int i = 0; i < n; ++i) {
        ipat[i] = itmp[i];
    }
}

} // namespace feff::path
