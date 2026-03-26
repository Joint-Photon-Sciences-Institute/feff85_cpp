#include "integration.hpp"
#include <cmath>

namespace feff::math {

void trap(const double x[], const double y[], int n, double& sum) {
    sum = y[0] * (x[1] - x[0]);
    for (int i = 1; i < n - 1; ++i) {
        sum += y[i] * (x[i + 1] - x[i - 1]);
    }
    sum += y[n - 1] * (x[n - 1] - x[n - 2]);
    sum /= 2.0;
}

void strap(const float x[], const float y[], int n, float& sum) {
    sum = y[0] * std::abs(x[1] - x[0]);
    for (int i = 1; i < n - 1; ++i) {
        sum += y[i] * std::abs(x[i + 1] - x[i - 1]);
    }
    sum += y[n - 1] * std::abs(x[n - 1] - x[n - 2]);
    sum /= 2.0f;
}

} // namespace feff::math
