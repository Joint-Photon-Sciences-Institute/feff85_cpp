#pragma once
// Loucks logarithmic radial grid utilities.
// Converted from: src/COMMON/xx.f
//
// The grid is x = log(r), with x(j) = -8.8 + (j-1)*0.05
// so r(j) = exp(x(j)). Index j is 1-based to match Fortran convention.

#include <cmath>

namespace feff::common {

/// Grid spacing in log-space
inline constexpr double grid_delta = 0.05;

/// Starting x value (x at j=1)
inline constexpr double grid_x0 = -8.8;

/// x-grid value at index j (1-based). x = log(r).
/// Replaces Fortran function xx(j).
inline double xx(int j) {
    return grid_x0 + (j - 1) * grid_delta;
}

/// Radial coordinate at grid index j (1-based). r = exp(x).
/// Replaces Fortran function rr(j).
inline double rr(int j) {
    return std::exp(xx(j));
}

/// Grid index of point immediately below radius r (1-based).
/// Replaces Fortran function ii(r).
inline int ii(double r) {
    return static_cast<int>((std::log(r) + (-grid_x0)) / grid_delta) + 1;
}

} // namespace feff::common
