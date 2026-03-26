#pragma once

// Distance calculations
// Converted from src/MATH/dist.f and src/MATH/sdist.f

#include <feff/types.hpp>

namespace feff::math {

// Distance between two 3D Cartesian points (double precision)
double dist(const double r0[3], const double r1[3]);

// Distance between two 3D Cartesian points (single precision)
float sdist(const float r0[3], const float r1[3]);

} // namespace feff::math
