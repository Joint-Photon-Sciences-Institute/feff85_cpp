#pragma once

// Trapezoidal integration
// Converted from src/MATH/trap.f and src/MATH/strap.f

namespace feff::math {

// Trapezoidal integration of y(x), result in sum. Double precision.
void trap(const double x[], const double y[], int n, double& sum);

// Trapezoidal integration, single precision, using abs(dx).
void strap(const float x[], const float y[], int n, float& sum);

} // namespace feff::math
