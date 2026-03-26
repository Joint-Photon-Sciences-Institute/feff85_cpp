// Complex energy grid construction for SCF-MT calculation.
// Converted from src/POT/grids.f

#include "grids.hpp"
#include <feff/constants.hpp>
#include <cmath>

namespace feff::pot {

void grids(double ecv, double xmu, int negx, int& neg,
           FeffComplex* emg, double* step, int nflrx)
{
    // eimmin = the lowest imaginary energy to search for Fermi level
    FeffComplex eimmin = coni * 0.05 / hart;

    int neg1 = (nflrx + 1) / 2;
    int neg3 = nflrx - 1;
    int neg2mx = negx - neg1 - neg3;

    // Never do calculations on real axis
    FeffComplex eim = eimmin;

    // Region 1: near ecv, imaginary part increases quadratically
    for (int i = 1; i <= neg1; ++i) {
        eim = eimmin * static_cast<double>(i * i);
        emg[i - 1] = ecv + eim;
    }
    step[nflrx - 1] = std::imag(eim) / 4.0;

    // Set energy step for integration, eim above real axis
    double de = std::imag(emg[neg1 - 1]) / 4.0;
    int neg2 = static_cast<int>(std::round((xmu - ecv) / de));
    if (neg2 > neg2mx) neg2 = neg2mx;
    if (neg2 < neg1) neg2 = neg1;
    de = (xmu - ecv) / neg2;

    // Region 2: horizontal path from ecv to xmu
    for (int i = neg1 + 1; i <= neg1 + neg2; ++i) {
        emg[i - 1] = emg[i - 2] + de;
    }

    neg = neg1 + neg2 + neg3;

    // Region 3: near xmu, imaginary part increases quadratically
    for (int i = 1; i <= neg3; ++i) {
        eim = eimmin * static_cast<double>((i + 1) * (i + 1)) / 4.0;
        if (i <= nflrx) step[i - 1] = std::imag(eim) / 4.0;
        emg[neg - i] = xmu + eim;
    }
}

} // namespace feff::pot
