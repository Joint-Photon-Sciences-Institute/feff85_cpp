#pragma once

// Physical and mathematical constants from FEFF
// Converted from src/HEADERS/const.h

#include <complex>

namespace feff {

// Mathematical constants
inline constexpr double pi = 3.1415926535897932384626433;
inline constexpr double one = 1.0;
inline constexpr double zero = 0.0;
inline constexpr double third = 1.0 / 3.0;
inline constexpr double raddeg = 180.0 / pi;

// Complex unit
inline const std::complex<double> coni(0.0, 1.0);

// kf = fa/rs, see Ashcroft and Mermin p.37
inline constexpr double fa = 1.919158292677512811;

// Physical constants
inline constexpr double bohr = 0.52917721067;     // Bohr radius in Angstrom
inline constexpr double ryd = 13.60569301;         // Rydberg in eV
inline constexpr double hart = 2.0 * ryd;          // Hartree in eV
inline constexpr double alpinv = 137.035999139;     // 1/alpha (fine structure)
inline constexpr double alphfs = 1.0 / alpinv;      // Fine structure constant

} // namespace feff
